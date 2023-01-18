function [Qlims,cIo,cOo,success] = conduitQfromRadius(C,Qbounds,varargin)
% val = conduitRadiusFromQ(C,lookupT,varargin)
% Given R and surface pressure condition (or general cI struct?), run 
% shooting search to get appropriate mass discharge rate (range) for
% conduit model (using Hajimirza conduit model, V6).
% IN:
%   C       = input conduit struct (ie from getConduitSource.m)
%   Qbounds = [2x1] bounding search values - valid R values lie within this
%             range, but Rbounds values themselves are not valid.
%
% Optional Name/Value pairs:
%   Qvalid      : scalar or 2x1, giving a single (or bounded) R values with
%                 valid conduit solutions (ie search in to out)
%   Zfailthresh : 10;   % Z (depth) threshold (m) (> Fails)
%   Mfailthresh : 0.95; % Mach number threshold (< Fails)
%   Pfailthresh : .1;  % Overpressure ratio threshold (> Fails)
% 
%   dQminScale  : 1e-3; % Tolerance for R search: Rtol = Ri*dRminscale
%   maxIter     : 10;      % Max iterations to narrow search
%   
%
% Nope:
%   lookupT = struct with 3d array fields Q, pf, a, Z, pm, M, cI0?
%             (fields corresponding to final (end) output vals of conduit model
%
% C Rowell, Jun 2021

Zfailthresh = 2;   % Z (depth) threshold in units of vent radius
Mfailthresh = 0.95; % Mach number threshold
Pfailthresh = .05;  % Fractional pressure difference threshold

dQminScale = 1e-3; % Tolerance for R search: Rtol = Ri*dRminscale
maxIter = 10;      % Max iterations to narrow search

% see also conduitQfromRadius.m?
% C Rowell, May 2021

% Tfile = '/Users/crrowell/Kahuna/data/glaciovolc/conduitSweeps/conduitV6_coarseSweep_n32_21-05-10.mat';
%% Parse input
%     C = getConduitSource;

    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    
%     addParameter(p,'searchVar','conduit_radius')
%     addParameter(p,'cI',C)
%     addParameter(p,'lookupT',Tfile)
%     addParameter(p,'Q',C.Q)
%     addParameter(p,'pf',C.pf)

    addParameter(p,'Qvalid',[])
    addParameter(p,'Zfailthresh',Zfailthresh)
    addParameter(p,'Mfailthresh',Mfailthresh)
    addParameter(p,'Pfailthresh',Pfailthresh)
    addParameter(p,'dQminScale',dQminScale)
    addParameter(p,'maxIter',maxIter)
    addParameter(p,'output',false)
    addParameter(p,'verbose',false)

    parse(p,varargin{:})
    par = p.Results;
       
    if par.verbose; tic; end
%     load(par.lookupT)

    % Make sure Qbounds have a minimum separation based on dQminScale 
    % (~4 iterations worth, or 1/4 maxIter)
    if diff(Qbounds) < (par.dQminScale)
        Qbounds = Qbounds + max(Qbounds).*par.dQminScale.*max([2.^round(par.maxIter/4) 2.^3]).*[-1 1];
    end
    
    %% Do the thing
    
    % Tests to bracket the success range
    % ==== Fine search - no successes, but bracketed ====
    if isempty(par.Qvalid)
        dQ = diff(Qbounds)/2;
        dQmin = max(Qbounds)*par.dQminScale;
        
%         passChecks = [false false];
        valid = false;
        
        iter = 0;
        Qlo = min(Qbounds);
        Qhi = max(Qbounds);
        if par.verbose; fprintf('Q search...\n  R: %.3f, Pf: %.3e, Q_lo: %.5e Q_hi: %.5e\n',C.conduit_radius,C.pf,Qlo,Qhi); end
        while and(dQ>dQmin,~valid) && iter<=par.maxIter*2
            iter=iter+1;
            C.Q = mean([Qlo Qhi]);
            cO = Conduit_flow_with_nucleation_V6(C);
            
            
            [Zcheck,UPcheck,Mcheck,Pcheck,flareCheck,valid] = checkConduitResult(cO,p.Zfailthresh,p,Mfailthresh,p.Pfailthresh);
                      
            if par.verbose
                fprintf('  --> I: %i, Q: %.5e, dQ: %.5e, Zc: %i, UPc: %i, Pc: %i, Mc: %i, Fc: %i, V: %i\n',...
                    iter,C.Q,dQ,Zcheck,UPcheck,Pcheck,Mcheck,flareCheck,valid)
            end
            
            % Q too low
            if (and(~Mcheck,Zcheck) && ~flareCheck)  % Z success, U fail, no flare
                Qlo = C.Q;
                dQ  = dQ/2;
            % Q too high
            elseif and(Mcheck,~Zcheck) ||... % U success, Z fail
                    and(and(~Mcheck,Zcheck), flareCheck) ||... % flare+Ufail+Zpass
                    (~Mcheck && ~Zcheck && Pcheck) || ...% Z fail, U fail, P pass (non-unique but should be too low in most useful cases)
                    and(~valid,~cO.Par.frag) % no Frag
                Qhi = C.Q;
                dQ  = dQ/2;
            end
            
        end
        if valid
            par.Qvalid = C.Q;
            cIo = C;
            cOo = cO;
        else
            warning('Could not converge on a valid Q with:\n\tdQ/dQmin=%.3f, nIter=%i, R=%.3f, pf=%.2f, [Qlo Qhi]=%.2f %.2f',...
                dQ/dQmin,iter,C.conduit_radius,C.pf/1e6,Qlo,Qhi)
            Qlims=[Qlo Qhi];
            cIo = C;
            if exist('cO','var')
                cOo = cO;
            else
                cIo.Q = mean(Qbounds);
                cOo = Conduit_flow_with_nucleation_V6(C);
            end
            success = false;
            if par.verbose; toc; end
            return
        end
    end
    
    % ==== Fine search - found successes ====
        for searchDir = [1 -1]  % --- Search up, then down ---
            if searchDir==1
                lastQsuccess = max(par.Qvalid); %ric(target1,jj,kk);
                lastQfail    = max(Qbounds);    %ric(iMax,jj,kk);
                passI        = 1;
            elseif searchDir==-1
                lastQsuccess = min(par.Qvalid); %ric(target2,jj,kk);
                lastQfail    = min(Qbounds);    %ric(iMin,jj,kk);
                passI        = 2;
            end

            dQ = abs(lastQsuccess-lastQfail)/2;
            dQmin = lastQsuccess*par.dQminScale;
            lastFlare = NaN;

            valid      = false;
            C.Q = lastQsuccess;
            iter = 0;
            if par.verbose; fprintf('Qbound...\n R: %.3f, Pf: %.3e, Qi: %.5f, dQi: %.5f\n',C.conduit_radius,C.pf,lastQsuccess,searchDir*dQ); end

            while or(dQ>dQmin,~valid) && iter<=par.maxIter
                iter = iter+1;

%                 lastRsuccess = C.conduit_radius;
                C.Q = C.Q + searchDir*dQ;
                cO = Conduit_flow_with_nucleation_V6(C);

                % Allow underpressure such that Pd + rho*v^2/2 ~ Pf
                
                [Zcheck,UPcheck,Mcheck,Pcheck,flareCheck,valid] = checkConduitResult(cO,p.Zfailthresh,p.Mfailthresh,p.Pfailthresh);
%                 Pcheck = and((cO.pm(end)/cO.Par.pf)<(1+par.Pfailthresh), ((cO.pm(end)+cO.U(end).^2*cO.rho_magma(end)/2)/cO.Par.pf)>(1-par.Pfailthresh));
%                 Mcheck     = and(cO.M(end)>par.Mfailthresh,cO.M(end)<1.1);
%                 passChecks(1) = or(Mcheck,... % Mach # Check (or P balance) (going up)
%                                 Pcheck); 
%                 passChecks(2) = cO.Z(end)<par.Zfailthresh; % Depth check (going down)
%                 flareCheck = cO.a(end)>C.conduit_radius;

                if par.verbose
                    fprintf('  --> I: %i, Q: %.5e, dQ: %.5e, Zc: %i, UPc: %i, Pc: %i, Mc: %i, Fc: %i, V: %i\n',...
                        iter,C.Q,searchDir*dQ,Zcheck,UPcheck,Pcheck,Mcheck,flareCheck,valid)
                end
                if valid 	 % Reduce step size and continue
                    dQ = dQ/2;
                    lastQsuccess = C.Q;
                    if par.output && searchDir==-1 % Generally want minimum Q for a non-flare solution
                        cIo = C;
                        cOo = cO;
                    end
                elseif flareCheck % Failed with conduit flare, so go back and reduce step size. SearchDir?
                    lastFlare = C.Q;
                    dQ = dQ/2;
                    C.Q = nanmin([lastQsuccess lastFlare]);
                    
                else % Go back if overshoot, reduce step size
                    dQ = dQ/2;
                    lastQfail = C.Q;
                    C.Q = nanmin([lastQsuccess lastFlare]);
                end

            end
            if searchDir==1
                Qlims(2) = lastQsuccess;
            elseif searchDir==-1
                Qlims(1) = lastQsuccess;
            end
            success = true;
        end
        if ~exist('cOo','var') && par.output
            cIo = C;
            cIo.Q = min(Qlims);
            cOo = Conduit_flow_with_nucleation_V6(cIo);
        end

        if par.verbose; fprintf('  OUT: Qf_lo: %.5e, Qf_hi: %.5e, t: %.4f s\n\n',Qlims(1),Qlims(2),toc);end
        if ~par.output; cIo=[]; cOo=[]; end
end