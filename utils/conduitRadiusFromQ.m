function [Rlims,cIo,cOo,validOutcomeCodes,allOutcomeCodes] = conduitRadiusFromQ(C,Rbounds,varargin)
% [Rlims,cIo,cOo,validOutcomeCodes,allOutcomeCodes] = conduitRadiusFromQ(C,Rbounds,varargin)
% Given Q and surface pressure condition (or general cI struct?), run 
% shooting search to get appropriate radius (range) for
% conduit model (using Hajimirza conduit model, V6).
% IN:
%   C       = input conduit struct (ie from getConduitSource.m)
%   Rbounds = [2x1] bounding search values - valid R values lie within this
%             range, but Rbounds values themselves are not valid.
%               TODO: Make this optional?
%
% Optional Name/Value pairs:
%   Rvalid      : scalar or 2x1, giving a single (or bounded) R values with
%                 valid conduit solutions (ie search in to out)
%
%   dRminScale  : 1e-3; % Tolerance for R search: Rtol = Ri*dRminscale
%
%   maxIter     : 10;      % Max iterations to narrow search
%
%   unboundedSearchMode : T/F - forces a step by step search through
%                          Rbounds parameter space. Rbounds are searched
%                          systematically rather than adaptively.
%                          The search occurs from center (mean) outwards
%                          and stops when a valid solution is found.
%                          Leaving Rbounds = [] forces this option, in
%                          which case Rbounds are guessed as 
%                          [0.6-1.4] .* extrapVentRadius(C.Q);
%
%   unboundedMaxIter    : Sets step size for unbounded search. 
%                         Default is 41, ie Rbounds = 0.6:0.02:1.4;
%
% C Rowell, May 2021

dRminScale = 1e-3; % Tolerance for R search: Rtol = Ri*dRminscale
maxIter = 10;      % Max iterations to narrow search

% see also conduitQfromRadius.m?
% C Rowell, May 2021
%    - UPDATED Dec 2023

%% Parse input
%     C = getConduitSource;

    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    
    addParameter(p,'Rvalid',[])
    addParameter(p,'dRminScale',dRminScale)
    addParameter(p,'maxIter',maxIter)
    addParameter(p,'output',false)
    addParameter(p,'unboundedSearchMode',false)
    addParameter(p,'unboundedMaxIter',41)
    addParameter(p,'verbose',false)

    parse(p,varargin{:})
    par = p.Results;
       
    if par.verbose; tic; end
%     load(par.lookupT)

    if isempty(Rbounds)
        par.unboundedSearchMode = true;
        Rbounds = extrapVentRadius(C.Q).*[0.6 1.4];
    end

    % Make sure Rbounds have a minimum separation based on dRminScale 
    % (~4 iterations worth, or 1/4 maxIter)
    if diff(Rbounds) < (par.dRminScale)
        Rbounds = Rbounds + max(Rbounds).*par.dRminScale.*max([2.^round(par.maxIter/4) 2.^3]).*[-1 1];
    end
    
    allOutcomeCodes = [];
    validOutcomeCodes = [0 0];
    %% Do the thing 1) Coarse search with no prior guess on Rbounds and no valid result
    
    if par.unboundedSearchMode && isempty(par.Rvalid)
        % Set up search range and search order
        Rlo = min(Rbounds);
        Rhi = max(Rbounds);
        searchScale = linspace(Rlo,Rhi,par.unboundedMaxIter);
        dR = mean(diff(searchScale));
        
        % Set index order to search outwards from mean, alternating
        % direction
        searchIndexOrder = zeros(par.unboundedMaxIter,1);
        searchIndexOrder(1:2:end) = ceil(par.unboundedMaxIter/2):-1:1;
        searchIndexOrder(2:2:end) = (ceil(par.unboundedMaxIter/2)+1):1:par.unboundedMaxIter;
        
        if par.verbose 
            fprintf('Unbounded search for initial valid result...guessing bounds...\n   Max %i tests: Q: %.2e, Pf: %.3e, R_lo: %.5f R_hi: %.5f\n',...
                par.unboundedMaxIter,C.Q,C.pf,Rlo,Rhi); 
        end

        iter = 1;
        valid = false;
%         largestR = searchScale(searchIndexOrder(iter));
%         smallestR = searchScale(searchIndexOrder(iter));
        while ~valid && iter <= par.unboundedMaxIter
            
            C.conduit_radius = searchScale(searchIndexOrder(iter));

            [cO,~] = conduitFlowRun(C);
            valid = cO.Outcome.Valid;
            if ~ismember(cO.Outcome.Code,allOutcomeCodes)
                allOutcomeCodes = [allOutcomeCodes cO.Outcome.Code];
            end
            
            if par.verbose
                reportString(cO.Outcome,C,iter,dR)
            end
            
            iter=iter+1;

        end
        if cO.Outcome.Valid
            par.Rvalid = C.conduit_radius;
            cIo = C;
            cOo = cO;
            
            % Reset Rbounds here based on new info...
            % ...Found downwards
            if searchScale(searchIndexOrder(iter)) > searchScale(searchIndexOrder(1))
                % Alternating search direction means 2 steps ago was last fail
                Rbounds = [searchScale(searchIndexOrder(iter-2)) max(Rbounds)];
                
            % ...Found upwards
            elseif searchScale(searchIndexOrder(iter)) < searchScale(searchIndexOrder(1))
                Rbounds = [min(Rbounds) searchScale(searchIndexOrder(iter-2)) ];
            
            % ... first guess - don't revise Rbounds in this case
%             elseif searchScale(searchIndexOrder(iter)) == searchScale(searchIndexOrder(1))
                
            end
            
        else
            Rlims=[Rlo Rhi];
            cIo = C;
            if exist('cO','var')
                cOo = cO;
            else
                cIo.conduit_radius = mean(Rbounds);
                cOo = conduitFlowRun(C);
            end
%             success = false;
            if par.verbose
                warning('No valid radius found with:\n\tdR/dRmin=%.3f, nIter=%i, Q=%.2e, pf=%.2f, [Rlo Rhi]=%.2f %.2f',...
                    dR/Rlo,par.unboundedMaxIter,C.Q,C.pf/1e6,Rlo,Rhi)
                toc; 
            end
            return
        end
    end
    
    %% Tests to bracket the success range
    % ==== Fine search - no successes, but bracketed ====
    if isempty(par.Rvalid) && ~par.unboundedSearchMode
        dR = diff(Rbounds)/2;
        dRmin = max(Rbounds)*par.dRminScale;
        
        valid = false;
        
        iter = 0;
        Rlo = min(Rbounds);
        Rhi = max(Rbounds);
        if par.verbose; fprintf('Searching for initial valid result...\n  Q: %.2e, Pf: %.3e, R_lo: %.5f R_hi: %.5f\n',C.Q,C.pf,Rlo,Rhi); end
        while and(dR>dRmin,~valid) && iter<=par.maxIter*2
            iter=iter+1;
            C.conduit_radius = mean([Rlo Rhi]);
            
            [cO,~] = conduitFlowRun(C);
            valid = cO.Outcome.Valid;
            if ~ismember(cO.Outcome.Code,allOutcomeCodes)
                allOutcomeCodes = [allOutcomeCodes cO.Outcome.Code];
            end
            
            if par.verbose
                reportString(cO.Outcome,C,iter,dR)
            end
            
            % Too high
            if and(~cO.Outcome.Choked,cO.Outcome.DepthFlag) && ~cO.Outcome.Flared
                Rhi = C.conduit_radius;
                dR  = dR/2;
                
            % Too low
            elseif and(cO.Outcome.Choked,~cO.Outcome.DepthFlag) ||... % U success, Z fail
                    and(and(~cO.Outcome.Choked,cO.Outcome.DepthFlag), cO.Outcome.Flared) ||... % flare+Ufail+Zpass
                    (~cO.Outcome.Choked && ~cO.Outcome.DepthFlag && cO.Outcome.PressureBalanced) || ...% Z fail, U fail, P pass (non-unique but should be too low in most useful cases)
                    and(~cO.Outcome.Valid,~cO.Outcome.Frag) % no Frag
                
                Rlo = C.conduit_radius;
                dR  = dR/2;
%             else
%                 error('Could not assess R adjustment condition.')
            end
            
        end
        if cO.Outcome.Valid
            par.Rvalid = C.conduit_radius;
            cIo = C;
            cOo = cO;
        else
            warning('Could not converge on a valid radius with:\n\tdR/dRmin=%.3f, nIter=%i, Q=%.2e, pf=%.2f, [Rlo Rhi]=%.2f %.2f',...
                dR/dRmin,iter,C.Q,C.pf/1e6,Rlo,Rhi)
            Rlims=[Rlo Rhi];
            cIo = C;
            if exist('cO','var')
                cOo = cO;
            else
                cIo.conduit_radius = mean(Rbounds);
                cOo = conduitFlowRun(C);
            end
%             success = false;
            if par.verbose; toc; end
            return
        end
    end
    
    %% ==== Fine search - found successes ====

        for searchDir = [1 -1]  % --- Search up, then down ---
            if searchDir==1
                lastRsuccess = max(par.Rvalid); %ric(target1,jj,kk);
                lastRfail    = max(Rbounds);    %ric(iMax,jj,kk);
                passI        = 1;
                searchStr    = 'upper';
            elseif searchDir==-1
                lastRsuccess = min(par.Rvalid); %ric(target2,jj,kk);
                lastRfail    = min(Rbounds);    %ric(iMin,jj,kk);
                passI        = 2;
                searchStr    = 'lower';
            end
%                 lastRfail = ric(iMax,jj,kk);
            dR = abs(lastRsuccess-lastRfail)/2;
            dRmin = lastRsuccess*par.dRminScale;
            lastFlare = NaN;

            valid      = false;
            C.conduit_radius = lastRsuccess;
            iter = 0;
            if par.verbose
                fprintf('Seeking R %s bound with limits: R_lo: %.5f, Rhi: %.5f\n  Q: %.2e, Pf: %.3e, Ri: %.5f, dRi: %.5f\n'...
                    ,searchStr,min(Rbounds),max(Rbounds),C.Q,C.pf,lastRsuccess,searchDir*dR); 
            end

            while or(dR>dRmin,~valid) && iter<=par.maxIter
                iter = iter+1;

                C.conduit_radius = C.conduit_radius + searchDir*dR;
                cO = conduitFlowRun(C);

                if par.verbose
                    reportString(cO.Outcome,C,iter,dR)
                end
                
                if cO.Outcome.Valid 	 % Reduce step size and continue
                    dR = dR/2;
                    validOutcomeCodes = updateValidOutcomeCodes(C.conduit_radius,...
                        cO.Outcome,validOutcomeCodes,lastRsuccess,searchDir);
                    lastRsuccess = C.conduit_radius;

                    if par.output && searchDir==1
                        cIo = C;
                        cOo = cO;
                        
                    end
                elseif cO.Outcome.Flared && cO.Outcome.DepthFlag % Failed, but conduit flares, so continue from here. DepthFlag req't added provisionally, Dec 2023. SearchDir? 
                    lastFlare = C.conduit_radius;
                    dR = dR/2;
                    
                else % Go back if overshoot, reduce step size
                    dR = dR/2;
                    lastRfail = C.conduit_radius;
                    C.conduit_radius = nanmax([lastRsuccess lastFlare]);
                end

            end
            if searchDir==1
                Rlims(2) = lastRsuccess;
%                 outcomeCode(2) = 
            elseif searchDir==-1
                Rlims(1) = lastRsuccess;
            end
%             success = true;
            if ~exist('cOo','var') && par.output
                cIo = C;
                cIo.conduit_radius = Rlims(2);
                cOo = conduitFlowRun(C);
            end
        end
        if par.verbose; fprintf('  OUT: Rf_lo: %.5f, Rf_hi: %.5f, t: %.4f s\n\n',Rlims(1),Rlims(2),toc);end
        if ~par.output; cIo=[]; cOo=[]; end
end

function reportString(Outcome,C,iter,dR)
    if Outcome.Valid; checkStr = '√'; elseif Outcome.Failed; checkStr = 'F'; else; checkStr = 'x';end
    fprintf(' %s -> I: %i, dR= %.5f, R= %.3f: %s\n',...
            checkStr,iter,dR,C.conduit_radius,Outcome.reportString)

end

function validOutcomeCodes = updateValidOutcomeCodes(thisR,Outcome,validOutcomeCodes,lastRsuccess,searchDir)
    if Outcome.Valid
        if searchDir == 1 && or(validOutcomeCodes(2) <= 0 , thisR > lastRsuccess)
            validOutcomeCodes(2) = Outcome.Code;
            
        elseif searchDir == -1 && or(validOutcomeCodes(1) <= 0 , thisR < lastRsuccess)
            validOutcomeCodes(1) = Outcome.Code;
            
        end
    end
end