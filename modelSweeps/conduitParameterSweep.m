function [dat,varMin,varMax,outcomeCodes] = conduitParameterSweep(cI_common,sweepParams,outDir,varargin) % runCoarse,runFine) % OR conduitPressureSweep?
% Conduit sweep function, conduit V7
% INPUT:
%   cI          : common input parameters for all model runs
%   sweepParams : Parameters to sweep over as a struct.
%                  Limit to 1-2 parameters max for now? maybe 3?.
%                  Struct fieldnames = parameter name
%                     substruct: .values  = [2x1]
%                     substruct: .n       = number of tests, this dimension
%   outDir      : Directory to save data. Will not save to disc if empty.
%
% OPTIONAL NAME/VALUE ARGS:
%   cores       : number of cores to run. Using >1 core will force 
%                 verbose = false within the conduit search function to
%                 avoid output chaos.
%   searchVar   : [OPTIONAL] For a given run:
%                   'R' = [Default] search over conduit radius, fixing mass flux
%                   'Q' = search over mass flux, fixing conduit radius
%                   ** ONLY 'R' is currenly implemented
%
%   descriptor  : text tag to add to save file names (so avoid spaces etc)
%
%   verbose :  T/[F], spit out all the details
%
% THESE OPTIONS AREN'T IN YET
%   runCoarse   : [T]/F - coarse solution search with fixed steps in R/Q
%                  (purpose is to find or bracket initial valid solutions)
%
%   runFine     : T/[F] - fine search with adaptive steps in R/Q
%                  (purpose is to get fine-tuned bounds on valid solutions)
%                   *** REQUIRES coarseFile if runCoarse=F? Or valid
%                   solutions input?
%                   % TODO: replace this with Rbounds values, scale, or function?
%
%   coarseFile  : path to save file for previously created output of a 
%                 coarse search. Automatically sets runCoarse = false
%
%   lookHarder  : Integer. DEFAULT = NaN. If no valid solutions are found in
%                  coarse search, refine the search and run another set 
%                  (will slow things down). 
%                   Typical value = 2, which means halve previous grid
%                   spacing and look in the previously untested values.
%                   3 = 1/3 spacing, etc.
%
%
%   PLUS additional params for conduitRadiusFromQ:
%       dRminScale, maxIter, unboundedMaxIter 
%
%  ---> Shooting runs to find the range of viable conduit solutions for
%  "dry" runs, aka Hw=0
%       --> Conduit stop conditions: Pm~Pf, U~0, M>~.999
%       --> Check for valid conduit outcomes
%       --> Require Z(end)~0 AND either (a) choked or (b) pm(end) = pf
%           --> So watch out for runs where stop condition is M>=.999
%  ---> Test a range of vent pressures (~water depths) on the base case to
%   assess response, adjusting MER or radius (phi0 or Z0 as options?)
%
% PROCEDURE (should write as accessible functions to refine full look-up table):
%     (1) "runCoarse": Coarse brute force sweep with a coarse step size to find
%           approximate bounding locations for conduit radii.
%           --> Generate a coarse lookup table in (params of interest) R,MER,Pf
%           --> CURRENTLY SAVES ALL CONDUIT OUTPUT, FILE CAN BE VERY LARGE
%               '-> Coarse search really only needs summary data (step 2)
%
%     (2) "getParams": Get summary values from coarse search, sufficient to
%           run fine search.
%           --> Could combine with (1) if you don't care to save coarse
%           output
%
%     (3) "runRefine": Use initial coarse table to:
%           (a) refine precision bounds on R for each MER,Pf 
%               --> (or go the other way and refine precision MER for given R,Pf)?
%
%   NOT YET FULLY IMPLEMENTED
%           (b) interpolate guesses for arbitrary R, MER, Pf from fine
%           search lookup table.
%
%
% C Rowell, December 2023
% Uses Hajimirza conduit model, updated to V7
disp('======================================')
disp('  Setting up conduit model sweep...')

    if nargin<2 || isempty(outDir)
        saverun = false;
    end
    
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;

    addParameter(p,'cores',         1)
%     addParameter(p,'fixed',         'R')
%     addParameter(p,'runCoarse',     true)
%     addParameter(p,'runFine',       false)
%     addParameter(p,'coarseFile',    '')
    addParameter(p,'descriptor',    '')
%     addParameter(p,'lookHarder',    false)
    addParameter(p,'verbose',       false)
    % PLUS conduitRadiusFromQ/conduitQFromRadius params...

    parse(p, varargin{:})
    runPars = p.Results;
    
    if ~exist(outDir,'dir')
        error('Output directory not found.')
    end
    
    if isstruct(sweepParams)
        [parNames,sweepVectors,sweepDims,nRuns] = parseSweepParamsStruct(sweepParams);
    end
    fprintf('\t%i Total Runs\n',nRuns)
    for pi=1:length(parNames)
        parStr(pi) = string(sprintf('%i (%s)',sweepDims(pi),parNames{pi}));
    end
    fprintf('\t[%s]\n',join(parStr,' x '))
    fprintf('\t%i Core(s)\n\n',runPars.cores)
    
    % Initialize output arrays
    disp('  Populating input parameters...')

    genericCI = getConduitSource();
    dat(nRuns) = struct('cI',genericCI,'cO', conduitFlowRun(genericCI));
    dat = reshape(dat,sweepDims);
    failMsg(nRuns) = struct('msg',[],'runi',[],'sweepPars',[]);
    failMsg = reshape(failMsg,sweepDims);
    modelFail  = false(sweepDims);

    % Populate all input params - doing this outside parfor for now, should
    % be reasonably fast.
    dat = populateConduitInputs(dat,cI_common,parNames,sweepVectors,sweepDims,nRuns);
    
    % Would use these to adapt R values from completed runs
%     completedValidRuns = false(dsz);
    varMin          = zeros(sweepDims);
    varMax          = zeros(sweepDims);
    varMinCodes     = zeros(sweepDims);
    varMaxCodes     = zeros(sweepDims);
    allCodes        = cell(sweepDims);
   
    if ~isempty(runPars.descriptor)
        rtitle = ['Conduit sweep: ' strrep(runPars.descriptor,'_',' ')];
    else
        rtitle = ['Conduit sweep'];
    end
    
    if and(runPars.cores>1, runPars.cores<=20)
        % parpool here
        multithread = true;
        parpool(runPars.cores);
        q = parallel.pool.DataQueue;
        afterEach(q, @nUpdateWaitbar);
        prog = 0;
    else
        multithread = false;
        prog = 0;
    end
    
    disp('Running sweep...')
    startTime = tic;
    h = waitbar(0, sprintf('%i / %i',prog,nRuns), 'Name', rtitle);
   
%     if ~runPars.verbose; textprogressbar('  Sweep: ');end
    if multithread
    
        parfor ii=1:nRuns

            tic
            % TODO: input optional to limit Rbounds or add Rvalid, other
            %       params for conduitRadiusFromQ
            try
                [Rlims,dat(ii).cI,dat(ii).cO,validCodes,allCodes{ii}] = conduitRadiusFromQ(dat(ii).cI,[],'verbose',false,'output',true);
            catch ME
                disp('da fuq')  % We shouldn't be erroring out at this point, as error catching should happen inside the search
                dat(ii).cO = ConduitOutcome.getErrorOutcomeFields;
                dat(ii).cO.Outcome = ConduitOutcome(ME);
                Rlims = [NaN NaN];
                validCodes = [NaN NaN];
                allCodes{ii} = dat(ii).cO.Outcome.Code;
            end
            varMin(ii) = Rlims(1);
            varMax(ii) = Rlims(2);
            varMinCodes(ii) = validCodes(1);
            varMaxCodes(ii) = validCodes(2);
            
%             rTime = toc; % For verbose reporting...
           

        % TODO: Adapting Rbounds here based on the previous run would
        % be helpful - get new average, tighten down the range, optimize
        % search spacing. Will work well since we are sweeping....
        %   - the clever move is to find adjacency in completed runs?
        %   - but hard to query in parfor loop...?
        %   - and tricky because we need to eval if adjacency is close
        %   enough...
            % UPDATE R values and use adjacency lookup?    

%             if cO.Outcome.Valid
%                 completedValidRuns(ii) = true;
%             end

            % TODO: finish adding in relevant code tidbits following
            % hydrovolcMC
            
            modelFail(ii) = dat(ii).cO.Outcome.Failed;
            if modelFail(ii)
                failMsg(ii).runi = ii;
                failMsg(ii).ME = dat(ii).cO.Outcome.Exception;
            end
            send(q, ii);
        end
        delete(gcp)

    else
        
        for ii=1:nRuns

            tic
            try
                [Rlims,dat(ii).cI,dat(ii).cO,validCodes,allCodes{ii}] = conduitRadiusFromQ(dat(ii).cI,[],'verbose',runPars.verbose,'output',true);
            catch ME
                dat(ii).cO = ConduitOutcome.getErrorOutcomeFields;
                dat(ii).cO.Outcome = ConduitOutcome(ME);
                Rlims = [NaN NaN];
                validCodes = [NaN NaN];
                allCodes{ii} = dat(ii).cO.Outcome.Code;
            end
            varMin(ii) = Rlims(1);
            varMax(ii) = Rlims(2);
            varMinCodes(ii) = validCodes(1);
            varMaxCodes(ii) = validCodes(2);
%             rTime = toc; % For verbose reporting
            
            modelFail(ii) = dat(ii).cO.Outcome.Failed;
            if modelFail(ii)
                failMsg(ii).runi = ii;
                failMsg(ii).ME = dat(ii).cO.Outcome.Exception;
            end
            waitbar(ii/nRuns, h, sprintf('%i / %i',ii,nRuns));
        end
    end
    close(h)
    disp('   ...Done!')
    toc(startTime)
    % TODO Get summary data... (optionally)
    outcomeCodes.minR       = varMinCodes;
    outcomeCodes.maxR       = varMaxCodes;
    outcomeCodes.all        = allCodes;
    outcomeCodes.modelFail  = modelFail;
    outcomeCodes.failMsg    = failMsg;
    
    
    % Save output
    if ~isempty(outDir)
        parStr = '';
        for pi=1:length(parNames)
            parStr = [parStr  sprintf('_%s_%i',parNames{pi},sweepDims(pi))];
        end
        
        if ~isempty(runPars.descriptor)
            saveStr = sprintf('%s_%s_%in%s',datestr(now,'yyyy-mm-dd'),runPars.descriptor,nRuns,parStr);
        else
            saveStr = sprintf('%s_%in%s',datestr(now,'yyyy-mm-dd'),nRuns,parStr);
        end
        fprintf('Saving output:\n\tDir: %s\n\tFile: %s\n',outDir,saveStr)
        save(fullfile(outDir,saveStr),'cI_common','dat','sweepParams','runPars','varMin','varMax','outcomeCodes')
        
    end
    disp('')
    
    function nUpdateWaitbar(~)
        prog = prog + 1;
        waitbar(prog/nRuns, h, sprintf('%i / %i',prog,nRuns));
%         p = p + 1;
    end
end

function dat = populateConduitInputs(dat,cI,parNames,sweepVectors,sweepDims,nRuns)
    for ni = 1:length(sweepDims)
        % Repmat vectors to get full list
        newshape = ones(1,length(sweepDims)); newshape(ni) = length(sweepVectors{ni});
        repvec = sweepDims; repvec(ni) = 1;
        sweepVectors{ni} = repmat(reshape(sweepVectors{ni},newshape),repvec);
    end
    
    for ni = 1:nRuns
        cIn = cI;
        for pi = 1:length(parNames)
            cIn.(parNames{pi}) = sweepVectors{pi}(ni);
        end

        dat(ni).cI = getConduitSource(cIn);
    end
end

function [parNames,sweepVectors,sweepDims,nRuns] = parseSweepParamsStruct(sweepParams)

    parNames = fieldnames(sweepParams);
    nPar = length(parNames);
    assert(nPar<=3,'Sweeping over more than 3 params simultaneously is not supported.')
    sweepDims = zeros(1,nPar);
    sweepVectors = cell(1,nPar);

    for ni = 1:nPar
        assert(ismember(parNames{ni},fieldnames(getConduitSource())),...
            sprintf('Sweep parameter name %s is not a recognized conduit input.',parNames{ni}))

        if isstruct(sweepParams.(parNames{ni}))
            sweepDims(ni) = sweepParams.(parNames{ni}).n;
            sweepVectors{ni} = linspace(sweepParams.(parNames{ni}).range(1),sweepParams.(parNames{ni}).range(2),sweepDims(ni))';

        elseif isvector(sweepParams.(parNames{ni}))
            sweepDims(ni) = length(sweepParams.(parNames{ni}));
            sweepVectors{ni} = sweepParams.(parNames{ni});
            if ~iscolumn(sweepVectors{ni})
                sweepVectors{ni} = sweepVectors{ni}';
            end

        else
            error('Input format of sweepParams field "%s" not recognized.',parNames{ni})
        end
    end
    
        nRuns = prod(sweepDims);

    if nRuns>5e5
        % Comment this out if you're sure...
        error('Current sweep params will perform %i runs. Be sure...',nRuns) 
    elseif nRuns>1e4
        error('Current sweep params will perform %i runs.',nRuns)
    end
end


