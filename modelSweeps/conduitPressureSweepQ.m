% Conduit pressure sweep V2: Fixed R and Pf, adjust MER
%   Find appropriate MER over Zw range given atmospheric profile and vent 
%   altitude, for a range of water depths, Zw
%       1) Get MER range for initial surface runs - find correct R for this base case
%           --> Use initial pressure sweep and adjust radius a bit as needed
%       2) Get MER adjustment for each base case as a function of Zw
%           --> Again interpolate from initial pressure sweep (fine sweep
%           as bounds, interpolated value as initial guess)
%
% C Rowell, June 2021
% Uses Hajimirza conduit model, V6
clear all; close all

% codeDir = '~/code/research-projects/glaciovolc/glaciovolc-dev/';
% dataDir = '/Users/crrowell/Kahuna/data/glaciovolc/conduitSweeps/';
codeDir = 'C:\Users\crowell\Documents\GitHub\glaciovolc\glaciovolc-dev';
dataDir = 'C:\Users\crowell\Kahuna\data\glaciovolc\conduitSweeps';

addpath(genpath(codeDir))

interpFile = fullfile(dataDir,'conduitV6_fineSweep_n64821_21-05-16_compressedV4_fineV2_rInterpV1.mat');

baseRunFile = fullfile(dataDir,'conduitV6_Zw0_baseRuns_Grimsvotn_V1.mat');
outFile     = 'conduitV6_%s_Q_from_R-vs-Zw_Grimsvotn_V1.mat';

% Setup params
Zw = [0:10:500]; % Water depths of interest
Rselect = 4:41; % Subset original MER/R sweep
Qminthresh    = 5e4; % Do not run (auto fail) any below this value

rho_l = 1000;
g     = 9.81;

% ---------- Atmo profiles and vent elevations ----------
% Tropical
% conIn.atmo  = fullfile(codeDir,'1Dplume_DB2012/atm_ERAreanalysis_Tungarahua2014_01_absWind.mat');
% conIn.vh0   = 0;

% High lat
conIn.atmo = fullfile(codeDir,'1Dplume_DB2012/atm_ERAreanalysis_Grimsvotn2011_01_absWind.mat');
conIn.vh0  = 1700;

% ----------- Search params -------------
% Conduit run failure thresholds
ZfailScale  = 2;   % Z (depth) threshold in conduit radii
Mfailthresh = 0.95; % Mach number threshold
Pfailthresh = .04;  % Over-/underpressure threshold

% Tolerances
dRminScale  = 0.2e-3;
dQminScale  = 0.2e-3;
maxIter     = 10;

% Workflow
getBaseRuns = false;
getQfromR   = true;
%%  Do the thing
% load(atmoF)
% atmo = atmprofile(~isnan(atmprofile(:,1)),:);
% 
% pfh = 10.^interp1(atmo(:,2),log10(atmo(:,1)*100),vh0,'pchip','extrap') + rho_l*g*Zw;

if getBaseRuns
    % (1) Select radii at surface pressure as fixed reference
    load(interpFile)
    cI0 = cI; clear cI
    cS0 = cS; clear cS
    
    Rsurf = Rmax(:,1);
%     pi = 10;
%     qi = 5;

    conIn.Zw = 0;
    % cI = getConduitSource(conIn);
    pf_new = zeros(length(Rsurf),1);

    for ii = 1:length(Rsurf)
        conIn.conduit_radius = Rsurf(ii);
        cIf = getConduitSource(conIn);
        pf_new(ii) = cIf.pf;

    end
    % Estimate radii from interpolation
    [P0,Q0] = meshgrid(pf,MER);
%     R_guess = griddata(P0(cS0.success),Q0(cS0.success),Rmax(cS0.success),pf_new,MER);
    rMxFun = scatteredInterpolant(P0(cS0.success),Q0(cS0.success),Rmax(cS0.success),'linear');
    rMnFun = scatteredInterpolant(P0(cS0.success),Q0(cS0.success),Rmin(cS0.success),'linear');
    R_guess = rMxFun(pf_new,MER);

    % Rbscale = max(Rsurf(2:end)./Rsurf(1:end-1));
    % Rbounds = R_guess.*[1/Rbscale Rbscale];
%     Rbounds = [cS0.Rmin(:,1) cS0.Rmax(:,1)] + diff([cS0.Rmin(:,1) cS0.Rmax(:,1)],[],2).*[-.5 .5];
    Rbounds = sort([rMnFun(cS0.Rmin(:,1),MER) R_guess] + diff([rMnFun(cS0.Rmin(:,1),MER) R_guess],[],2).*[-.5 .5],2);
    
    % Get refined radius bounds
    cIf = getConduitSource(conIn);
    Rmax_base = R_guess*0;
    Rmin_base = R_guess*0;
    success_base   = R_guess*0;

    for ii = 1:length(Rsurf)
        cIf.conduit_radius = R_guess(ii);
        cIf.Q = MER(ii);

        Zfailthresh = ZfailScale*cIf.conduit_radius;

        cO = Conduit_flow_with_nucleation_V6(cIf);
        [~,~,~,~,~,valid] = checkConduitResult(cO,Zfailthresh,Mfailthresh,Pfailthresh);
        if valid
            Rvalid = cIf.conduit_radius;
        else
            Rvalid = [];
        end    

        [Rrange,cIi,cO,success_base(ii)] = conduitRadiusFromQ(cIf,Rbounds(ii,:),'Rvalid',Rvalid,'maxIter',maxIter,...
            'Zfailthresh',Zfailthresh,'Mfailthresh',Mfailthresh,'Pfailthresh',Pfailthresh,...
            'dRminScale',dRminScale,'verbose',false,'output',true);

        cI_base(ii) = table2struct(cIi);
        Rmax_base(ii) = max(Rrange);
        Rmin_base(ii) = min(Rrange);

        fprintf(' %i/%i ... Q: %.2e, Pf: %.2f, Rb_lo: %.3f\tRb_hi: %.3f\tRf_lo: %.4f\tRf_hi: %.4f\tV?: %i\n',...
        ii,length(Rsurf),cIi.Q,cI_base(ii).pf/1e6,Rbounds(ii,1),Rbounds(ii,2),min(Rrange),max(Rrange),success_base(ii))

    end
    save(baseRunFile,'cI_base','Rmax_base','Rmin_base','success_base')
else
    load(baseRunFile)
end

%% (2) Get Q as function of pf(Zw)

if getQfromR
    load(interpFile)
    load(baseRunFile)
    
    [P0,Q0] = meshgrid(pf,MER);
    Rvalid = ~isnan(Rmax);
    qFmin = scatteredInterpolant(P0(Rvalid),Rmax(Rvalid),Q0(Rvalid),'linear','nearest');
    qFmax = scatteredInterpolant(P0(Rvalid),Rmin(Rvalid),Q0(Rvalid),'linear','nearest');
    
    % Select subset
    Rmax = Rmax(Rselect,:);
    Rmin = Rmin(Rselect,:);
    Rmax_base = Rmax_base(Rselect);
    Rmin_base = Rmin_base(Rselect);
    success_base = success_base(Rselect);
    MER  = MER(Rselect);
    cI = cI(Rselect,:);
    cI_base = cI_base(Rselect);
    Rvalid = ~isnan(Rmax);
    
%     cS0  = cS0(Rselect,:);
%     cI0  = cI0(Rselect,:);
    
%     Qmax = zeros(length(MER),length(Zw));
%     Qmin = Qmax;
%     Qsuccess = false(size(Qmin));
    cIq(length(MER),length(Zw)) = cI_base(end);
    
    dd = NaN(length(MER),length(Zw));
    cQ = struct('R',Rmax_base,'Zw',Zw,'pf',dd,'Qmin',dd,'Qmax',dd,'Ztop',dd,'M',dd,'pm',dd,'rho_b',dd,...
        'U',dd,'amax',dd,'dPdt_fr',dd,'phi_g',dd,'Z_fr',dd,'U_fr',dd,...
        'pfr',dd,'success',false(size(dd)),'qSrchBounds',repmat(dd,[1 1 2]),...
        'zFailThresh',dd); %,'checks',zeros(length(MER),6));


    
    ct = 0;
    nSearches = numel(dd);
    for ii = 1:length(MER)
        for jj = 1:length(Zw)
            ct  = ct+1;
            cIf = cI_base(ii);
            cIf = rmfield(cIf,'pf');
            cIf.Zw = Zw(jj);
            cIf = getConduitSource(cIf);
            
%             Qbounds = MER([qi-1 qi+1]);
%             Qguess  = interp1(Rmax(:,pi),MER,Rmax_base(ii),'spline');
            
%             Qguess  = griddata(P0(Rvalid),Rmax(Rvalid),Q0(Rvalid),cIf.pf,cIf.conduit_radius);
%             Qguess_max = griddata(P0(Rvalid),Rmin(Rvalid),Q0(Rvalid),cIf.pf,cIf.conduit_radius);
            
%             qFmin = scatteredInterpolant(P0(Rvalid),Rmax(Rvalid),Q0(Rvalid),'linear');
%             qFmax = scatteredInterpolant(P0(Rvalid),Rmin(Rvalid),Q0(Rvalid),'linear');
            Qguess = qFmin(cIf.pf,cIf.conduit_radius);
            Qguess_max = qFmax(cIf.pf,cIf.conduit_radius);
            
            cIf.Q   = Qguess;

            Zfailthresh = ZfailScale*cIf.conduit_radius;

            if Rvalid(ii,jj) && Qguess>=Qminthresh
                cO = Conduit_flow_with_nucleation_V6(cIf);
                [~,~,~,~,~,valid] = checkConduitResult(cO,Zfailthresh,Mfailthresh,Pfailthresh);
                if valid
                    Qvalid = cIf.Q;
                else
                    Qvalid = [];
                end

                Qbounds = sort([Qguess Qguess_max]);
                Qbounds = Qbounds + diff(Qbounds).*[-.5 .5];

                [Qrange,cIi,cO,success] = conduitQfromRadius(cIf,Qbounds,'Qvalid',Qvalid,'maxIter',maxIter,...
                    'Zfailthresh',Zfailthresh,'Mfailthresh',Mfailthresh,'Pfailthresh',Pfailthresh,...
                    'dQminScale',dQminScale,'verbose',false,'output',true);

                cIq(ii,jj) = table2struct(cIi);            
    %         cQ.checks(jj,:) = [Zpass,UPpass,Mpass,Ppass,Fcheck,valid];
    %             cQ.R(ii)    = cIi.conduit_radius;
                cQ.Qmin(ii,jj) = min(Qrange);
                cQ.Qmax(ii,jj) = max(Qrange);
                cQ.success(ii,jj) = success;
                cQ.Ztop(ii,jj)    = cO.Z(end);
                cQ.M(ii,jj)       = cO.M(end);
                cQ.pf(ii,jj)      = cO.Par.pf;
                cQ.pm(ii,jj)      = cO.pm(end);
                cQ.rho_b(ii,jj)   = cO.rho_magma(end);
                cQ.U(ii,jj)       = cO.U(end);
                cQ.amax(ii,jj)    = cO.a(end);
                cQ.phi_g(ii,jj) = cO.porosity(end);
                cQ.Z_fr(ii,jj)  = cO.Par.Zf;          

                if cO.Par.frag
                    fi                  = find(cO.Z==cO.Par.Zf,1,'first');
                    cQ.dPdt_fr(ii,jj)   = cO.dPdt(fi); % Peak dPdt at fragmentation
                    cQ.pfr(ii,jj)       = cO.pm(fi);
                    cQ.U_fr(ii,jj)      = cO.U(fi);     % Velocity at fragmentation
                else
                    cQ.dPdt_fr(ii,jj)   = max(cO.dPdt); % Peak dPdt at fragmentation
                    cQ.pfr(ii,jj)       = NaN; %cO.pm(fi);
                    cQ.U_fr(ii,jj)      = NaN; %cO.U(fi);     % Velocity at fragmentation
                end
            else
                Qbounds = [NaN NaN];
            end
            
            fprintf(' %i/%i ... R: %.4f,  Pf: %.2f,  Qb_lo: %.4e,  Qb_hi: %.4e,  Qf_lo: %.4e,  Qf_hi: %.4e,  V?: %i\n',...
                ct,nSearches,cQ.R(ii),cQ.pf(ii,jj)/1e6,Qbounds(1),Qbounds(2),cQ.Qmin(ii,jj),cQ.Qmax(ii,jj),cQ.success(ii,jj))
            
%             if jj==2
%                 disp('Pause here')
%             end
        end
    end
    outName = fullfile(dataDir,sprintf(outFile,datestr(now,'YYYY-mm-dd')));
    fprintf('Saving:\n\t%s\n',outName)
    save(outName,'cI_base','Rmax_base','Rmin_base','success_base','cIq','cQ',...
        'Rselect','Zfailthresh','Mfailthresh','Pfailthresh')
    
end

%%

% for ii = qi %1:length(R)
% 
% %     dd = zeros(size(R,1),1); %dd = zeros(length(R),length(pf));
%     
%     cIf = cI(ii,pi);
%     cIf.conduit_radius = R(ii);
%     % Test Q search
%     Qbounds = MER([qi-1 qi+1]);
%     Qguess  = interp1(Rmax(:,pi),MER,R(ii),'spline');
%     cIf.Q   = Qguess;
%     Zfailthresh = ZfailScale*cIf.conduit_radius;
% 
%     cO = Conduit_flow_with_nucleation_V6(cIf);
%     [~,~,~,~,~,valid] = checkConduitResult(cO,Zfailthresh,Mfailthresh,Pfailthresh);
%     if valid
%         Qvalid = cO.Q;
%     else
%         Qvalid = [];
%     end
%     [Qrange,cIi,cO,success] = conduitQfromRadius(cIf,Qbounds,'Qvalid',Qvalid,'maxIter',maxIter,...
%         'Zfailthresh',Zfailthresh,'Mfailthresh',Mfailthresh,'Pfailthresh',Pfailthresh,...
%         'dQminScale',dQminScale,'verbose',true,'output',true);
%     
%     
%     % Run a test sweep
%     [P,Q] = meshgrid(pf,MER);
% %     IFun = scatteredInterpolant(P(cS.success),Rmax(cS.success),Q(cS.success));
%     Qr = griddata(P(cS.success),Rmax(cS.success),Q(cS.success),pf,R(ii)*ones(size(pf)));
%     Qvec    = linspace(min(Qr),max(Qr),length(pf));
%     dd      = zeros(size(Qvec)); 
%     cQ = struct('Q',Qvec,'Ztop',dd,'M',dd,'pm',dd,'rho_b',dd,...
%         'U',dd,'amax',dd,'dPdt_fr',dd,'phi_g',dd,'Z_fr',dd,'U_fr',dd,...
%         'pfr',dd,'success',false(size(dd)),'qSrchBounds',repmat(dd,[1 1 2]),...
%         'zFailThresh',dd,'checks',zeros(length(Qvec),6));
%     
%     cQ.Ztop = dd;
%     cQ.M    = dd;
%     cQ.Q    = dd;
%     cQ.pf   = dd;
%     cQ.pm   = dd;
%     cQ.amax = dd;
%     for jj = 1:length(Qvec)
%         cIf.Q = Qvec(jj);
%         cO = Conduit_flow_with_nucleation_V6(cIf);
%         Zthresh = ZfailScale*cIf.conduit_radius;
%         [Zpass,UPpass,Mpass,Ppass,Fcheck,valid] = checkConduitResult(cO,Zthresh,Mfailthresh,Pfailthresh);
%         
%         cQ.checks(jj,:) = [Zpass,UPpass,Mpass,Ppass,Fcheck,valid];
%         cQ.R           = R(ii);
%         cQ.Ztop(jj)    = cO.Z(end);
%         cQ.M(jj)       = cO.M(end);
%         cQ.pf(jj)      = cO.Par.pf;
%         cQ.pm(jj)      = cO.pm(end);
%         cQ.rho_b(jj)   = cO.rho_magma(end);
%         cQ.U(jj)       = cO.U(end);
%         cQ.amax(jj)    = cO.a(end);
% 
% 
%         
%     end
% end

%%
% pi = 1;

% figure
% ax(1) = subplot(4,1,1);
% plot(Qvec,cQ.Ztop./(cQ.R*ZfailScale),'.-')
% hold on
% yl = ylim;
% plot([1;1]*Qrange,yl'*[1 1],'r')
% ylabel('Z/Z_{thresh}')
% axis tight
% 
% ax(2) = subplot(4,1,2);
% plot(Qvec,cQ.M,'.-');
% hold on
% yl = ylim;
% plot([1;1]*Qrange,yl'*[1 1],'r')
% ylabel('M')
% axis tight
% 
% ax(3) = subplot(4,1,3);
% plot(Qvec,cQ.pm./cQ.pf,'.-')
% hold on
% yl = ylim;
% plot([1;1]*Qrange,yl'*[1 1],'r')
% ylabel('K')
% axis tight
% 
% ax(4) = subplot(4,1,4);
% plot(Qvec,cQ.amax./cQ.R)
% hold on
% yl = ylim;
% plot([1;1]*Qrange,yl'*[1 1],'r')
% ylabel('a_v/a_0')
% axis tight
% 
% linkaxes(ax,'x')