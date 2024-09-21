function [dat,varMin,varMax,outcomeCodes] = refineConduitSweep(datfile,outDir,varTol,pool,qcplot)
% [dat,varMin,varMax,outcomeCodes] = refineConduitSweep(datfile,outDir,varTol,pool,qcplot)
% this function will take in a conduit parameter sweep, and use the search
% variable (R or Q) for successful runs to interpolate a refined search
% range for the failed runs. In many cases this should allow successful
% runs where before there was none.
%   
% PARAMETERS:
%   datfile = path to saved sweep output file
%   outDir  = new path for output - required to save output to disk
%   varTol  = relative tolerance to apply to interpolated input bounds
%       -> default is 0.1
%       -> search bounds for each run will be:
%     [varMax_interpolated * (1-varTol)  varMax_interpolated * (1+varTol)]
%   pool = optional parallel pool object

if nargin<2
    outDir = [];
end
if nargin<3 || isempty(varTol)
    varTol = 0.1;
end
if nargin<4
    make_parpool = true;
else
    make_parpool = false;
end
if nargin<5 || isempty(qcplot)
    qcplot = false;
end
default_n_cores = 4;

load(datfile,'dat','outcomeCodes','cI_common','runPars','sweepParams',...
    'varMax','varMin');
% dat = D.dat;
% cI_common = D.cI_common;
% outcomeCodes = D.outcome


%% Get booleans of all valid, invalid, and failed results

    [modelValid,modelInvalid,~] = checkValid(outcomeCodes);
    
    validMinVar = varMin(modelValid);
    validMaxVar = varMax(modelValid);

%% Get X,Y vectors for sweeps
    fn = fieldnames(sweepParams);
    assert(length(fn)==2,'Refine function is currently only built for 2D sweeps')

    x = linspace(sweepParams.(fn{1}).range(1),sweepParams.(fn{1}).range(2),sweepParams.(fn{1}).n);
    y = linspace(sweepParams.(fn{2}).range(1),sweepParams.(fn{2}).range(2),sweepParams.(fn{2}).n);


    % Swap X and Y if we need to
    if and(size(varMin,1)==length(x), size(varMin,2)==length(y)) && length(x)~=length(y)
        xt = x;
        x = y;
        y = xt;
    elseif length(x)==length(y)
        warning('X and Y vector lengths are ambiguous.')
    else
        assert(and(size(varMin,1)==size(X,1), size(varMin,2)==size(X,2)), 'Sweep variable X and Y sizes do not match outputs.')
    end
    
    [X,Y] = meshgrid(x,y);
    xi = X(modelValid);
    yi = Y(modelValid);
    
    xInv = X(modelInvalid);
    yInv = Y(modelInvalid);

%% Get interpolants
    intMethod = 'natural';
    extMethod = 'nearest';
    Fmin = scatteredInterpolant(xi,yi,validMinVar,intMethod,extMethod);
    Fmax = scatteredInterpolant(xi,yi,validMaxVar,intMethod,extMethod);
    varMinEst = Fmin(X,Y) .* (1-varTol);
    varMaxEst = Fmax(X,Y) .* (1+varTol);

%% Run the refined search for all invalid, non-failed results...

% nRuns = sum(modelInvalid(:));
nRuns = numel(dat);
sweepDims = size(dat);
invSimNum = find(modelInvalid);

% We run this in parallel by default for now
    if make_parpool
        pool = parpool(default_n_cores);
    end
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateWaitbar);
    prog = 0;

    disp('Running refined sweep search...')
    startTime = tic;
    if ~isempty(runPars.descriptor)
        rtitle = ['Conduit refined sweep: ' strrep(runPars.descriptor,'_',' ')];
    else
        rtitle = ['Conduit refined sweep'];
    end
    h = waitbar(0, sprintf('%i / %i',prog,nRuns), 'Name', rtitle);

    % Initializing arrays
    dat_new(nRuns)  = struct('cI',dat(end).cI,'cO', dat(end).cO);
    dat_new         = reshape(dat_new,sweepDims);
    failMsg(nRuns)  = struct('msg',[],'runi',[],'sweepPars',[],'ME',[]);
    failMsg         = reshape(failMsg,sweepDims);
    modelFail       = false(sweepDims);
    varMinCodes     = zeros(sweepDims);
    varMaxCodes     = zeros(sweepDims);
    allCodes        = cell(sweepDims);
    
    parfor ii=1:nRuns
        if modelInvalid(ii)
            [Rlims,dat_new(ii).cI,dat_new(ii).cO,validCodes,allCodes{ii}] = ...
                conduitRadiusFromQ(dat(ii).cI,[varMinEst(ii) varMaxEst(ii)],...
                'verbose',false,'output',true);
%             [Rlims,cI,cO,validCodes,aCodes] = ...
%                 conduitRadiusFromQ(dat(ii).cI,[],...
%                 'verbose',false,'output',true);

            % Vars that need updating: dat, outcomeCodes, varMin, varMax
%             dat(ii).cI = cI;
%             dat(ii).cO = cO;
            varMin(ii) = Rlims(1);
            varMax(ii) = Rlims(2);   
            varMinCodes(ii) = validCodes(1);
            varMaxCodes(ii) = validCodes(2);
%             allCodes{ii} = aCodes;
            
            modelFail(ii) = dat(ii).cO.Outcome.Failed;
            if modelFail(ii)
                failMsg(ii).runi = ii;
                failMsg(ii).ME = dat(ii).cO.Outcome.Exception;
            end            
        else
            dat_new(ii) = dat(ii);

        end
        send(q, ii);
    end
    
    close(h)
    disp('   ...Done!')
    fprintf(' -> ')
    toc(startTime)
    
    % Overwrite with new outcome code results
    outcomeCodes.minR(invSimNum)       = varMinCodes(invSimNum);
    outcomeCodes.maxR(invSimNum)       = varMaxCodes(invSimNum);
    outcomeCodes.all(invSimNum)        = allCodes(invSimNum);
    outcomeCodes.modelFail(invSimNum)  = modelFail(invSimNum);
    outcomeCodes.failMsg(invSimNum)    = failMsg(invSimNum);
    
    dat = dat_new;
    
    [newModelValid,newModelInvalid,~] = checkValid(outcomeCodes);

    fprintf(' -> Refined search results:\n\t Invalid: %i --> %i\n\t Valid:   %i --> %i\n',...
        sum(modelInvalid(:)),sum(newModelInvalid(:)),sum(modelValid(:)),sum(newModelValid(:)) ) 
    %% QC plot of interpolated bounds
    if qcplot
        figure
        % subplot(1,2,1)
        scatter3(xi,yi,validMinVar,'o','filled')
        hold on
        scatter3(xi,yi,validMaxVar,'o','filled')
        scatter3(X(:),Y(:),varMinEst(:),'o')
        % title('Search variable minima')
        % xlabel('X')
        % ylabel('Y')

        % subplot(1,2,2)
        hold on
        scatter3(X(:),Y(:),varMaxEst(:),'o')
        title('Refined search variable minima/maxima')
        xlabel('X')
        ylabel('Y')
        legend({'Known minima','Known maxima','New search minima','New search maxima'})


        % Get plottable codes
        plotCodes = getCodeSummaryArray(outcomeCodes); % Get corrected code array
        [cmap,cax,cticks,clabels,outcomeIndex,~] = outcomeColorMap(plotCodes,true, false);

        % Make the plot
        figure
        imagesc(x,y, outcomeIndex )
        colormap(gca,cmap)
        set(gca,'YDir','normal')
        caxis(cax)
        xlabel(fn{1})
        ylabel(fn{2})
        ctlcb = colorbar(gca,'location','eastoutside');
        caxis(cax);
        ctlcb.Ticks = cticks;
        ctlcb.TickLabels = clabels;

    end

%% Writing output
if ~isempty(outDir)
    [iDir, fname, ext] = fileparts(datfile);
    fprintf(' -> Saving refined output file to: \n\t %s\n',[fullfile(outDir,fname) ,ext])
    save(fullfile(outDir,fname),'cI_common','dat','sweepParams','runPars','varMin','varMax','outcomeCodes')

end

%% Waitbar fun
    function nUpdateWaitbar(~)
        prog = prog + 1;
        waitbar(prog/nRuns, h, sprintf('%i / %i',prog,nRuns));
%         p = p + 1;
    end
end

function [modelValid,modelInvalid,modelFailed] = checkValid(outcomeCodes)
% Check sweep outcomeCodes struct to get boolean arrays of valid and
% invalid conduit runs.
    modelFailed = outcomeCodes.modelFail;

    nSuccess = cellfun(@(x) sum(x>0),outcomeCodes.all); % Check for more than 1 success type
    if any(nSuccess(:) > 1)
        warning('Multiple successful solutions found:\n\t -> %s', datfile)
        pause  %Hold up
    end
    modelValid = nSuccess>0;
    modelInvalid = and(nSuccess<=0, ~modelFailed);
end