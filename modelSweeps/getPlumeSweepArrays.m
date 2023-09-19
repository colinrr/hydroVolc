function [qA, thresh] = getPlumeSweepArrays(dat,MCflag)
% Retrieves summary data for a set of 1DPlume model (only) runs. 
%   Note: for simplicity, this script assumes some parameters do not change
%   between runs inside the same file
%       > (e.g.  conduit proxy flag, atmo profile, particle size range)
%   dat = structure array from parameter sweep or monte carlo run
%   MC  : true = monte carlo run, false = param sweep
%     saveOut = [true/false] - will save output file to the same directory
%                           as input
%             = 'char'  - saveOut is used as the output save directory
%
%  OUTPUT:
%     qA     : structure of key parameter summary arrays
%     thresh : behavior thresholds. Param sweep only, not currently
%               implemented for plume-only runs
%

if nargin<2
    MCflag = false;
end

sigma_t = 1e6; % Melt tensile strength for strain rate frag criterion estimate
sLen = 10; % Smooth length for threshold curves
sLen2 = 5;

% Default Conduit run failure thresholds
ZfailScale  = 2;   % Z (depth) threshold in conduit radii
Mfailthresh = 0.95; % Mach number threshold
Pfailthresh = .05;  % Over-/underpressure threshold
%%

% A single sweep should use the same atmo profile for all realizations for now

% Detects first non-empty run
ei = 1;
while and(isempty(dat(ei).pI),ei<length(dat))
    ei = ei+1;
end


% proxy = dat(ei).cI.proxy;

% Fields to retrieve. "_ex" are extra fields to calculate
% if ~proxy
%     cIf = {'phi_frag','conduit_radius','Zw','pf','Q','vh0','T'};
%     cOf = {'pm','Z','a','Cm','rho_magma','U','M','C','cin'}; cOf_ex = {'K','frag','choke','dPdt_fr','u_fr','mu_fr','pm_fr','pg_fr','Zfr','cin_fr','drdt_fr','r_fr'};% K, frag params (dPdt,u,mu,Z,etc), Par flags (frag, choke, Zf)
% else
%     cIf = {'phi_frag','conduit_radius','Zw','pf','Q','vh0','T','n_ec','n_0'};
%     cOf = {'pm','pg','Z','a','Cm','rho_magma','rho_g','U','M','C','porosity','n','K'}; cOf_ex = {'K','frag','choke','Q','Q0','rho_melt'};% K, frag params (dPdt,u,mu,Z,etc), Par flags (frag, choke, Zf)
% end

% dOf = {'Ld','a','u','rho_d','c','pd','chi'}; dOf_ex = {'SSA','SAv'}; % SSA, SAv
% wOf = {'Lj','failedPlume','m','m_s','m_v','m_l','alpha','xv','z','T'}; wOf_ex = {'T0'};
pIf = {'n_0','xv_0','r_0','rho_B0','rho_g0','rho_m','T0','T_g','u_0','vh0'}; pIf_ex = {'SSA','SAv','useDecompressLength','rho_a0','Ri0'}; % SSA, SAv
pOf = {'hm','hb','collapse','rm','m','m_d','m_l','m_v','rho_B','theta','m_0','u','angle'}; pOf_ex = {'m_s','SSA_hm','SSA_hb','SAv_hb','nsi','nsi_tp','m_w_tp','m_s_tp','m_hb','m_l_hb','m_v_hb','m_l_tp','m_v_tp','m_tp','u_hb','rho_hb'}; % SSA_m, SAv_m

%%  QC Calcs

% Check success levels: totally blank, conduit fail, MWI fail, plume fail


ff = fieldnames(dat); 
qA.QCexist = zeros(length(dat),length(ff));
qA.QClevel = zeros(length(dat),1);
for ii=1:length(dat(:))
    
    % 1) check for completely empty stuctures
    for jj=1:length(ff)
        qA.QCexist(ii,jj) = ~isempty(dat(ii).(ff{jj}));
    end
    
%     % 2) Conduit fail?
% %     [~,ci] = ismember('cI',ff);
%     if qA.QCexist(ii,1) && isfield(dat(ii).cI,'proxy') && dat(ii).cI.proxy
%         qA.QClevel(ii,1) = 1;
%     elseif qA.QCexist(ii,1)
%         [~,~,~,~,~,qA.QClevel(ii,1)] = checkConduitResult(dat(ii).cO,dat(ii).cI.conduit_radius*ZfailScale,Mfailthresh,Pfailthresh);
%     end
    
    % 3) Plume input fail?
%     if qA.QCexist(ii,1) && ~dat(ii).wO.failedPlume
%         qA.QClevel(ii,1) = 1;
%     end
    
    % 3) Plume output fail?
    if qA.QCexist(ii,2) && isfield(dat(ii).pO,'u') && (dat(ii).pO.u(end).*sin(dat(ii).pO.angle(end)))<1
        qA.QClevel(ii,1) = 1;
    end
        
end



%%  Initialize arrays

% fp = zeros(size(dat));
% xv_end = fp;
% zw = fp;

qSz = size(dat);

% Prepopulate arrays
% for ff=1:length(cIf)
%     qA.cI.(cIf{ff}) = NaN(qSz);
% end
% for ff=1:length(cOf)
%     qA.cO.(cOf{ff}) = NaN(qSz);
% end
% for ff=1:length(cOf_ex)
%     qA.cO.(cOf_ex{ff}) = NaN(qSz);
% end
% for ff=1:length(dOf)
%     qA.dO.(dOf{ff}) = NaN(qSz);
% end
% for ff=1:length(dOf_ex)
%     qA.dO.(dOf_ex{ff}) = NaN(qSz);
% end
% for ff=1:length(wOf)
%     qA.wO.(wOf{ff}) = NaN(qSz);
% end
% for ff=1:length(wOf_ex)
%     qA.wO.(wOf_ex{ff}) = NaN(qSz);
% end
for ff=1:length(pIf)
    qA.pI.(pIf{ff}) = NaN(qSz);
end
for ff=1:length(pIf_ex)
    qA.pI.(pIf_ex{ff}) = NaN(qSz);
end
for ff=1:length(pOf)
    qA.pO.(pOf{ff}) = NaN(qSz);
end
for ff=1:length(pOf_ex)
    qA.pO.(pOf_ex{ff}) = NaN(qSz);
end

pi = find(qA.QClevel(:,1),1,'first');
qA.pO.nsi    = NaN([qSz size(dat(pi).pO.nsi,2)]);
qA.pO.nsi_tp = NaN([qSz size(dat(pi).pO.nsi,2)]);

%% LOOP AND RETRIEVE FIELDS
qA.atmo = dat(ei).pI.atmo;
qA.Rgsd = dat(ei).pI.Rgsd;
qA.atmo_fields = {'Pressure','Altitude ASL','Temperature','Wind_abs','RelHumidity'};
qA.htropo = findTPheight(qA.atmo(:,2)/1e3,qA.atmo(:,3))*1e3;

% 
for ii=1:numel(dat)
    [jj,kk] = ind2sub(qSz,ii);
    
    % --------------------- CONDUIT INPUT -------------------------
%     for ff=1:length(cIf)
%         if qA.QCexist(ii,1)  %~isempty(dat(ii).cI)
%             qA.cI.(cIf{ff})(ii) = dat(ii).cI.(cIf{ff})(end);
%         end
%     end
    
    % --------------------- CONDUIT OUTPUT ------------------------
%     if qA.QCexist(ii,2) %~isempty(dat(ii).cO)
%         for ff=1:length(cOf)
%             qA.cO.(cOf{ff})(ii) = dat(ii).cO.(cOf{ff})(end);
%         end
%         qA.cO.K(ii)       = qA.cO.pm(ii)./qA.cI.pf(ii);
%         qA.cO.frag(ii)    = dat(ii).cO.Par.frag;
%         qA.cO.choke(ii)   = dat(ii).cO.Par.choke;
%         
%         if ~proxy
%             qA.cO.Zfr(ii)     = dat(ii).cO.Par.Zf;
%             fi                = find(dat(ii).cO.Z==dat(ii).cO.Par.Zf,1,'first');
%             qA.cO.dPdt_fr(ii) = dat(ii).cO.dPdt(fi);
%             qA.cO.u_fr(ii)    = dat(ii).cO.U(fi);
%             qA.cO.mu_fr(ii)   = dat(ii).cO.mu(fi);
%             qA.cO.pm_fr(ii)   = dat(ii).cO.pm(fi);
%             qA.cO.pg_fr(ii)   = dat(ii).cO.pg(fi);
%             qA.cO.cin_fr(ii)  = dat(ii).cO.cin(fi);
%             qA.cO.drdt_fr(ii) = dat(ii).cO.drdt(fi);
%             qA.cO.r_fr(ii)    = dat(ii).cO.M1(fi)./dat(ii).cO.M0(fi);
%         else
%             qA.cO.Zfr(ii)     = NaN;
%             qA.cO.dPdt_fr(ii) = NaN; %dat(ii).cO.dPdt(fi);
%             qA.cO.u_fr(ii)    = NaN; %dat(ii).cO.U(fi);
%             qA.cO.mu_fr(ii)   = NaN; %dat(ii).cO.mu(fi);
%             qA.cO.pm_fr(ii)   = NaN; %dat(ii).cO.pm(fi);
%             qA.cO.pg_fr(ii)   = NaN; %dat(ii).cO.pg(fi);
%             qA.cO.cin_fr(ii)  = NaN; %dat(ii).cO.cin(fi);
%             qA.cO.drdt_fr(ii) = NaN; %dat(ii).cO.drdt(fi);
%             qA.cO.r_fr(ii)    = NaN; %dat(ii).cO.M1(fi)./dat(ii).cO.M0(fi);
%             qA.cO.Q(ii)       = dat(ii).cO.Par.Q;
%             qA.cO.Q0(ii)      = dat(ii).cO.Par.Q0;
%         end
%     end
%     
    
    % ---------------------- DECOMPRESSION OUTPUT ----------------------
%     if qA.QCexist(ii,4) %~isempty(dat(ii).dO)
%         for ff=1:length(dOf)
%             qA.dO.(dOf{ff})(ii) = dat(ii).dO.(dOf{ff})(end);
%         end
%         qA.dO.SSA(ii)  = sum(dat(ii).dO.nsi.*(3*dat(ii).pI.saScale./(dat(ii).dO.rhoi.*dat(ii).dO.Rgsd)));
%         qA.dO.SAv(ii) = qA.dO.SSA(ii) .* (1-dat(ii).dO.n_0).*dat(ii).cI.Q ./ (dat(ii).dO.u .* pi.* dat(ii).dO.a.^2);
%     end
    
    % ---------------------- MWI OUTPUT ------------------------
%     for ff=1:length(wOf)
%         if qA.QCexist(ii,4) && isfield(dat(ii).wO,wOf{ff})
%             qA.wO.(wOf{ff})(ii) = dat(ii).wO.(wOf{ff})(end);
%         end
%         if qA.QCexist(ii,4) && isfield(dat(ii).wO,'T')
%             qA.wO.T0(ii) = dat(ii).wO.T(1);
%         end
%     end
    
    % ---------------------- PLUME INPUT ------------------------
    if qA.QCexist(ii,1)
        qA.pI.useDecompressLength(ii) = dat(ii).pI.useDecompressLength;
%         if ~dat(ii).wO.failedPlume
        for ff=1:length(pIf)
            qA.pI.(pIf{ff})(ii) = dat(ii).pI.(pIf{ff})(end);
        end
        qA.pI.SSA(ii)  = sum(dat(ii).pI.nsi.*(3*dat(ii).pI.saScale./(dat(ii).pI.rhoi.*dat(ii).pI.Rgsd)));
        qA.pI.SAv(ii)  = qA.pI.SSA(ii) .* dat(ii).pO.m_s(1) ./ (dat(ii).pI.u_0 .* pi.* dat(ii).pI.r_0.^2);  
        qA.pI.rho_a0(ii) = dat(ii).pO.atmo.rho(1);
        [~,qA.pI.Ri0(ii)]  = getRichardsonProfile(dat(ii),[],'simple');
    end
    
    % ---------------------- PLUME OUTPUT ------------------------
    if qA.QCexist(ii,2) && ~isnan(dat(ii).pO.m_0) %isfield(dat(ii).pO,pOf{ff})
        for ff=1:length(pOf)
            qA.pO.(pOf{ff})(ii) = dat(ii).pO.(pOf{ff})(end);
        end
        
        qA.pO.SSA_hm(ii) = sum(dat(ii).pO.nsi(end,:)'.*(3*dat(ii).pI.saScale./(dat(ii).pI.rhoi.*dat(ii).pI.Rgsd)));
        if ~qA.pO.collapse(ii) % Bouyancy level props
            [~,zi] = min(abs(dat(ii).pO.z - qA.pO.hb(ii))); % Find LNB index
            [~,tpi] = min(abs((dat(ii).pO.z + dat(ii).pI.vh0) - qA.htropo));
            qA.pO.SSA_hb(ii)  = sum(dat(ii).pO.nsi(zi,:)'.*(3*dat(ii).pI.saScale./(dat(ii).pI.rhoi.*dat(ii).pI.Rgsd)));
            qA.pO.SAv_hb(ii) = qA.pO.SSA_hb(ii) .* dat(ii).pO.m_s(zi) ./ (dat(ii).pO.u(zi) .* pi.* dat(ii).pO.r(zi).^2);
            qA.pO.m_hb(ii)   = dat(ii).pO.m(zi);
            qA.pO.m_l_hb(ii) = dat(ii).pO.m_l(zi);
            qA.pO.m_v_hb(ii) = dat(ii).pO.m_v(zi);
            qA.pO.u_hb(ii)   = dat(ii).pO.u(zi);
            qA.pO.rho_hb(ii) = dat(ii).pO.rho_B(zi);
                
            % Tropopause props
            if (dat(ii).pO.z(end) + dat(ii).pI.vh0)>=qA.htropo % crossing tropopause?
                qA.pO.nsi_tp(jj,kk,:) = dat(ii).pO.nsi(tpi,:);
                qA.pO.m_s_tp(ii) = dat(ii).pO.m_s(tpi);
                qA.pO.m_l_tp(ii) = dat(ii).pO.m_l(tpi);
                qA.pO.m_v_tp(ii) = dat(ii).pO.m_v(tpi);
                qA.pO.m_tp(ii)   = dat(ii).pO.m(tpi);
                qA.pO.m_w_tp(ii) = dat(ii).pO.m_l(tpi) + dat(ii).pO.m_v(tpi);
            else
                qA.pO.m_s_tp(ii) = 0;
                qA.pO.m_l_tp(ii) = 0;
                qA.pO.m_v_tp(ii) = 0;
                qA.pO.m_tp(ii)   = 0;
                qA.pO.m_w_tp(ii) = 0;
            end
        end
        qA.pO.m_s(ii) = sum(dat(ii).pO.m_si(end,:));
        qA.pO.nsi(jj,kk,:) = dat(ii).pO.nsi(end,:);
    end
end
qA.pO.collapse(isnan(qA.pO.collapse)) = false;
%%  Calc regime threshold values

if ~MCflag  % Thresholds mostly relevant for parameter sweeps, which are not currently implemented for Plume-only runs
    qA.Q0 = qA.cI.Q(:,1); % NEEDS RECALCULATING FROM PLUME
    % logQ0 = log10(Q0);
%     qA.Zw = qA.cI.Zw(end,:);

    nQ  = length(qA.Q0);
%     nZw = length(qA.Zw);

%     LdZ_thresh    = NaN(nQ,1);
%     Ld0_thresh    = LdZ_thresh;
%     fPlume_thresh = LdZ_thresh;
%     LdLj_thresh   = LdZ_thresh;
%     dPdt_crit     = NaN(nQ,2);
%     M_thresh      = LdZ_thresh;
%     K_thresh      = LdZ_thresh;
    Clps_lo       = zeros(nQ,1);
    Clps_hi       = zeros(nQ,1);
    tp_Max        = NaN(nQ,1);
    tp_Min        = NaN(nQ,1);
    for ii=1:nQ
%         li = find((qA.cO.M(ii,:)>=0.95),1,'last'); 
%         if ~isempty(li); M_thresh(ii) = qA.cI.Zw(ii,li); end
% 
%         li = find((qA.cO.K(ii,:)<=1.05),1,'first'); 
%         if ~isempty(li); K_thresh(ii) = qA.cI.Zw(ii,li); end
% 
%         li = find(qA.dO.Ld(ii,:)<=qA.cI.Zw(ii,:),1,'first');
%         if ~isempty(li); LdZ_thresh(ii) = qA.cI.Zw(ii,li); end
% 
%         li = find(qA.dO.Ld(ii,:)<=0,1,'first');
%         if ~isempty(li); Ld0_thresh(ii) = qA.cI.Zw(ii,li); end
% 
%         li = find(qA.wO.failedPlume(ii,:),1,'first');
%         if ~isempty(li); fPlume_thresh(ii) =  qA.cI.Zw(ii,li); end
% 
%         li = find((qA.dO.Ld(ii,:)+qA.wO.Lj(ii,:))<=qA.cI.Zw(ii,:),1,'first'); 
%         if ~isempty(li); LdLj_thresh(ii) = qA.cI.Zw(ii,li); end

        li = find(~qA.pO.collapse(ii,:),1,'first')-1; 
        if ~or(isempty(li),li<1); Clps_lo(ii) = qA.cI.Zw(ii,li); end

        clps_chk = ~qA.pO.collapse & ~qA.wO.failedPlume;
        li = find(clps_chk(ii,:),1,'last')+1; 
        if ~(or(isempty(li),li>length(qA.Zw))); Clps_hi(ii) = qA.cI.Zw(ii,li); end

        % Tropopause injection?
        tpiMax = find((qA.pO.hb(ii,:) + qA.cI.vh0(ii,:) + qA.cI.Zw(ii,:))>=qA.htropo,1,'last');
        tpiMin = find((qA.pO.hb(ii,:) + qA.cI.vh0(ii,:) + qA.cI.Zw(ii,:))>=qA.htropo,1,'first');
        if ~isempty(tpiMax)
            tp_Max(ii) = qA.Zw(tpiMax);
            tp_Min(ii) = qA.Zw(tpiMin);
        end
    end
%     qA.cO.sigma_t = sigma_t;
% ls = qA.cO.dPdt_fr.*qA.cI.conduit_radius./qA.cO.u_fr.*(qA.cI.phi_frag./sigma_t);


%     qA.thresh.LdZ_thresh = smooth(LdZ_thresh,sLen);
%     qA.thresh.Ld0_thresh = smooth(Ld0_thresh,sLen);
%     qA.thresh.fPlume_thresh = smooth(fPlume_thresh,sLen);
%     qA.thresh.LdLj_thresh = smooth(LdLj_thresh,sLen);
%     qA.thresh.M_thresh = smooth(M_thresh,sLen);
    % Clps_lo = smooth(Clps_lo,sLength);
    qA.thresh.Clps_hi = smooth(Clps_hi,sLen);
    qA.thresh.Clps_lo = Clps_lo; %smooth(Clps_lo,sLength);
    qA.thresh.Clps_lo(Clps_lo==0) = NaN;

    qA.thresh.tp_Max = tp_Max;
    qA.thresh.tp_Min = tp_Min;
    qA.thresh.tp_Max(~isnan(tp_Max)) = smooth(tp_Max(~isnan(tp_Max)),sLen2);

%     thresh.LdZ      = LdZ_thresh;
%     thresh.Ld0      = Ld0_thresh;
%     thresh.fPlume   = fPlume_thresh;
%     thresh.LdLx     = LdLj_thresh;
%     thresh.M        = M_thresh;
    thresh.Clps_hi  = Clps_hi;
    thresh.Clps_lo  = Clps_lo;
    thresh.tp_Max   = tp_Max;
    thresh.tp_Min   = tp_Min;
else
    thresh = [];
end


%% 
% if saveOut
%     if ~exist('oDir','var')
%         [oDir,oName,oext] = fileparts(sweepFile);
%     else
%         [~,oName,oext] = fileparts(sweepFile);
%     end
%     oFile = fullfile(oDir,[oFileName oName oext]);
%     fprintf('Saving output summary:\n\t%s\n',oFile)
%     save(oFile,'qA','thresh')
% else
%     oFile = '';
% end

end


