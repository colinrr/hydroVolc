function varargout = plumeMC(pI0,MC,N,cores,randPars)
%   [summ,randPars,dat,qA,thresh] = hydroVolcMC(cI,pI,MC,n,cores)
%   Run Monte Carlo simulations for the plume model only.
% INPUT:
% % % % %    cI = fixed conduit input params
%    pI = fixed MWI/plume input params
%    MC = struct of randomization params.
%         Required fields: cI, pI (params to randomize for each of the
%                                    two main input structs)
%         Each parameter needs 2 sub-fields: 
%         1) MC.<"pI","cI">.<var_name>.dist : Distribution type w/ options:
%               'uniform', 'discrete','normal','function'
%         2) MC.<"pI","cI">.<var_name>.range : Range of values, depending
%               on "dist"
%                   'uniform'   = [1 x 2] = [low high]
%                   'discrete'  = vector of all possible values
%                   'normal'     = [1 x 2] = [mean std.dev]
%                   'function'  = text function to evaluate. When an input should
%                       be a function of another randomized input, pass text to 
%                       run 'eval' cmd, Eg. using cI,pI params:
%                            MC.cI.vh0.dist = 'function';
%                            MC.cI.vh0.range = '-cI.Zw';          
%                         (function handles weren't really working here)
%   
%           Unique cI fields: 'logQ'  : log10(Q), to get Q distribution
%                                       based on the log value.
%                             'a_var' : Specify range as a fraction of 
%                                        the optimum conduit choking radius. 
%                                        REQUIRES either cI.Q or cI.logQ
%
%   cores = number of cores to run
%   randPars = pre-filled table of random parameters. This is used to re-run
%               a full set of sims using the same parameters. 
%                  - In this case, cI0, pI0, and MC should all be the same from original run
%                  - N will be overwritten to size(randPars,1). 
%                  - Also check "atmo" file path in this case. 
%                  - Default is [] for a new run. 
%
%   output setting? full vs abbreviated model output?
%
% OUTPUT: 
%   summ     = table summary of run results
%   randPars = outputs a table of the randomized variables
%   dat      = [n x 1] struct containing full monte carlo simulation output.
% 
% If 4th output below is specified, summary arrays are
% automatically calculated:
%   qA       = summary output arrays of key params from getPlumeSweepArrays.m)
%
%
% NOTES
%   For paralellization, current set to max 20 cores. Adjust if needed.
%   C Rowell, Feb 2022

if nargin<5
    randPars = [];
end
if nargin<4
    cores = [];
end

% if isempty(cI0)
%     cI0 = struct();
% end
if isempty(pI0)
    pI0 = struct();
end

% if ~isfield(MC,'cI')
%     MC.cI = struct();
% end
if ~isfield(MC,'pI')
    MC.pI = struct();
end

% Conduit run failure thresholds
% ZfailScale  = 2;   % Z (depth) threshold in conduit radii
% Mfailthresh = 0.95; % Mach number threshold
% Pfailthresh = .05;  % Over-/underpressure threshold

%% Setup parpool?
if isempty(cores)
    cores = 1;
end
if and(cores>1, cores<=20)
    % parpool here
  multithread = true;
  parpool(cores);
else
    multithread = false;
end

% Pre-build output summary and data structs
% cfn        = fieldnames(MC.cI);
pfn        = fieldnames(MC.pI);

if isempty(randPars)
    clear randPars
    nRandPars  = length(pfn);

    rfn = [pfn];
    aa = ismember(rfn,'logQ');
    if any(aa); rfn{aa} = 'Q'; end
    aa = ismember(rfn,'a_var');
    if any(aa); rfn{aa} = 'conduit_radius'; end
    sargs = [rfn'; cell(1,nRandPars)]; 

    randPars(N)  = struct(sargs{:});
    reRunMC = false;
    
else
    % Overwrite inputs using pre-existing randPars table
    nRandPars = size(randPars,2);
    N         = size(randPars,1);
    randPars  = table2struct(randPars);
    reRunMC   = true;
    
end

% conCheck   = zeros(N,1);
modelFail  = zeros(N,1);
failMsg(N) = struct('ME',[],'runi',[],'randPars',[]);
dat(N,1)     = struct('pI',[],'pO',[]);
fn         = fieldnames(dat);

% rnames    = cell(1,nRandPars);
% rnames(:) = {'double'};
% randPars  = table('Size',[N nRandPars],'VariableTypes',rnames,'VariableNames',[cfn;pfn]);


% if isfield(cI0,'proxy')
%     if cI0.proxy
%         conType = 'proxy';
%     else
%         conType = 'full';
%     end
% else
%     conType = 'full';
% end
fprintf('Running 1DPLUME model Monte Carlo sweep:\n\tN = %i, %i Random Params, %i Core(s)\n',N,nRandPars,cores)

if multithread
    
    parfor ii=1:N
%         esplode
        tic
        if reRunMC
            [pI] = getMCreRunPars(pI0,MC,randPars(ii));
        else
            [pI,randPars(ii)] = populateParams(pI0,MC);
        end
        
        try
            pI = getPlumeSource(pI);
            [dat(ii).pO,dat(ii).pI] = hmodel(pI,true);
            modelFail(ii) = false;
            rTime = toc;
            reportOutputPar(dat(ii),N,ii,rTime,full(modelFail(ii)),randPars(ii))
            
        catch ME
            modelFail(ii) = true;
            failMsg(ii).runi = ii;
            failMsg(ii).randPars = randPars(ii);
            failMsg(ii).ME = ME;
            for fs = 1:length(fn)
                if ~isfield(dat(ii),fn{fs})
                    
                    switch fn{fs}
                        case 'pI'
                            dat(ii).pI = pI;
%                         case 'cI'
%                             dat(ii).cI = cI;
                        otherwise
                            dat(ii).(fn{fs}) = [];
                    end
                end
            end
        end
        
%         dat(ii) = dati;
        
        
%         randPars(ii) = rP;
        
%         dat(ii) = dati;
%         randPars(ii,:) = randParsi;
    end
    delete(gcp)
        
else
    
    ct = 0;
    for ii=1:N
        if reRunMC
            [pI] = getMCreRunPars(pI0,MC,randPars(ii)); % FIX F'N
        else
            [pI,randPars(ii)] = populateParams(pI0,MC); % FIX F'N
        end
        
        ct = ct+1;
        tic
        try
            pI = getPlumeSource(pI);
            [dat(ii).pO,dat(ii).pI] = hmodel(pI,true);
%             dat(ii) = runCoupledModel(cI,pI);
            modelFail(ii) = false;
            
            rTime = toc;
            reportOutput(dat(ii),N,ct,rTime,full(modelFail(ii)),randPars(ii))
            
        catch ME
            modelFail(ii) = true;
            failMsg(sum(modelFail)).runi = ii;
            failMsg(sum(modelFail)).randPars = randPars(ii);
            failMsg(sum(modelFail)).ME = ME;
            for fs = 1:length(fn)
                if ~isfield(dat(ii),fn{fs})
                    
                    switch fn{fs}
                        case 'pI'
                            dat(ii).pI = pI;
                        case 'cI'  % Can prob comment out
                            dat(ii).cI = cI;
                        otherwise
                            dat(ii).(fn{fs}) = [];
                    end
                end
            end
        end
        
               
%         dat(ii) = dati;
%         randPars(ii,:) = randParsi;
    end
    
end

%% Sort output?
    summ.modelFail = modelFail;
    summ.failMsg = failMsg;
    randPars = struct2table(randPars,'AsArray',true);

%% Automatically output summary arrays?
varargout{1} = summ;
varargout{2} = randPars;
varargout{3} = dat;
if nargout>=4
    [varargout{4}, ~] = getPlumeSweepArrays(dat,true);
end
% if nargout==5
%     varargout{5} = thresh;
% end

end

% Function to generate randomized model input
function [pI,randPars] = populateParams(pI,MC)
% Populate input structs using randomization parameters in MC struct
% randPars = vector of randomized variable results


    
    sn = {'pI'};
%     ff{1} = fieldnames(MC.cI);
    ff{1} = fieldnames(MC.pI);
    
%     fo{1} = cI;
    fo{1} = pI;
    
    randPars = struct();
    fungroups = cell(1,2); %Input function group
    
    for ss = 1:length(ff)
        if ~isempty(ff{ss})
            for ii = 1:length(ff{ss})
                fn = ff{ss}{ii};

    %             switch fn 
    %                 case 'logQ_rng'
    %                     fo{1}.Q = rand.*range(MC.logQ_rng) + min(MC.logQ_rng);
    %                     
    %                 case 'a_var'
    %                     a_peak = extrapVentRadius(cI.Q);
    %                     fo{1}.conduit_radius = a_peak.*(MC.a_var.*(2.*rand-1) + 1);
    %                     
    %                 otherwise

                assert(all(ismember({'dist','range'},fieldnames(MC.(sn{ss}).(fn)))),...
                    'Monte Carlo input: Variable "%s" is missing a required field',fn)
                dist = MC.(sn{ss}).(fn).dist;
                par  = MC.(sn{ss}).(fn).range;


                switch dist
                    case 'uniform'
                        fo{ss}.(fn) = range(par).*rand + min(par);

                    case 'discrete'
                        fo{ss}.(fn) = par(randi(length(par)));

                    case 'normal'
                        fo{ss}.(fn) = randn.*par(2) + par(1);

                    case 'function'
                        varfun.var = fn;
                        varfun.fun = par;
                        fungroups{ss} = [fungroups{ss} varfun];
                        fo{ss}.(fn) = [];
%                         fungroups{ss} = [fungroups{ss} {par}];
                        
                    otherwise
%                         if isa(par,'function_handle')
%                             fungroups{ss}.fun = par;
%                             fungroups{ss}.var = dist;
%                             tvar = eval(dist);
%                             fo{ss}.(fn) = range(tvar);
%                         else
                            warning('Distribution type "%s" not recognized, variable "%s" will not be randomized.',MC.(sn{ss}).(fn).dist,fn)
                            fo{ss}.(fn) = mean(par);
%                         end
                end        

                randPars.(fn) = fo{ss}.(fn);

            end
        end
    end

    % Special mass flux and conduit radii params % FIX THESE - ESP Q
    % Mass flux
%     if isfield(fo{1},'logQ')
%         fo{1}.Q = 10^fo{1}.logQ;
%         fo{1} = rmfield(fo{1},'logQ');
%         
%         randPars.Q = fo{1}.Q;
%         randPars = rmfield(randPars,'logQ');
%     end
%     % Conduit radius lookup
%     if isfield(fo{1},'a_var')
%         a_peak = extrapVentRadius(fo{1}.Q);
%         fo{1}.conduit_radius = a_peak.*fo{1}.a_var;
%         fo{1} = rmfield(fo{1},'a_var');
%         
%         randPars.conduit_radius = fo{1}.conduit_radius;
%         randPars = rmfield(randPars,'a_var');
% 
%     end
    
%     cI = fo{1};
    pI = fo{1};
    
    % Run any input params that are a function of random params
%     if ~isempty(fungroups{1})
%         for jj=1:length(fungroups{1})
%             cI.(fungroups{1}(jj).var) = eval(fungroups{1}(jj).fun);
% %             cI.(fungroups{1}(jj).var) = varfun();
%         end
%     end
    if ~isempty(fungroups{1})
        for jj=1:length(fungroups{1})
            pI.(fungroups{1}(jj).var) = eval(fungroups{1}(jj).fun);
        end
    end
%     randPars = struct2table(randPars);
    


end

function [pI] = getMCreRunPars(pI,MC,randPars)
% Populate input structs using previous MC run
% randPars is a single scalar struct here

%     cff = fieldnames(MC.cI);
    pff = fieldnames(MC.pI);

%     for ci = 1:length(cff)
% 
%         if strcmp(cff{ci},'logQ') % FIX
%     %         cI.Q = randPars(idx,'Q');
%             cff{ci} = 'Q';
%         elseif strcmp(cff{ci},'a_var')
%     %         cI.conduit_radius = randPars(idx,'conduit_radius');
%             cff{ci} = 'conduit_radius';
%         end
%         cI.(cff{ci}) = randPars.(cff{ci});
%     %     end
%     end

    for pi = 1:length(pff)
        pI.(pff{pi}) = randPars.(pff{pi});
    end
end

function reportOutput(dat,N,ct,rTime,modelFail,randPars)
%         if ~isfield(dat.cI,'Zw')
%             Zw = 0;
%         else
%             Zw = dat.cI.Zw;
%         end
%         if ~isempty(dat.cO)
%             cF = dat.cO.Par.conValid;
%             if ~isnan(cF)
%                 cF = ~cF;
%             end
%         else
%             cF = NaN;
%         end
%         if ~isempty(dat.wO)
%             fP = dat.wO.failedPlume;
%         else
%             fP = NaN;
%         end
        if ~isempty(dat.pO)
            clp = dat.pO.collapse;
            if isfield(dat.pO,'hm') && ~isnan(dat.pO.hm); hm  = dat.pO.hm/dat.pO.atmo.ztropo;else; hm=NaN; end
        else
            clp = NaN;
            hm = NaN;
        end


        % Spit out some output %% ADD back Q reporting
        fprintf('  %i/%i:   runtime = %.0f s, modFail = %i, Clps = %i, Z/Ztp = %.2f\n',...
            ct,N,rTime,modelFail,clp,hm)
        
        if modelFail
            disp('Failed run:')
            disp(randPars)
        end

end

function reportOutputPar(dat,N,ct,rTime,modelFail,randPars)
%         if ~isfield(dat.cI,'Zw')
%             Zw = 0;
%         else
%             Zw = dat.cI.Zw;
%         end
%         if ~isempty(dat.cO)
%             cF = dat.cO.Par.conValid;
%             if ~isnan(cF)
%                 cF = ~cF;
%             end
%         else
%             cF = NaN;
%         end
%         if ~isempty(dat.wO)
%             fP = dat.wO.failedPlume;
%         else
%             fP = NaN;
%         end
        if ~isempty(dat.pO)
            clp = dat.pO.collapse;
            if isfield(dat.pO,'hm') && ~isnan(dat.pO.hm); hm  = dat.pO.hm/dat.pO.atmo.ztropo;else; hm=NaN; end
        else
            clp = NaN;
            hm = NaN;
        end


        % Spit out some output
        fprintf('  %i/%i:   runtime = %.0f s, modFail = %i, Clps = %i, Z/Ztp = %.2f\n',...
            ct,N,rTime,modelFail,clp,hm)
        
        if modelFail
            disp('Failed run:')
            disp(randPars)
        end

end