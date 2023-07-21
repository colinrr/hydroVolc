function [pI2,wO] = MWIv2(cI,cO,pI)
% [pI2,wO] = MWIv2(cI,cO,pI,dO)
%   Run near-vent magma-water interaction model.
% IN: 
%       cI    conduit input param struct
%       cO    conduit output struct
%       pI    plume input param struct
%
% OUT:  pI2   updated plume source parameter struct for plume model
%       wO    output struct for magma-water interaction model
%             FIELDS: 
%
%
% C Rowell, Jun 2021


%% SETUP
%     global xv
    P.xv_crit     = 0.04; % Dryness fraction below which to stop run
    
    % Properties and Constants
    P.g         = 9.81;
    P.Tr        = 273.15+1; % Could also be Tw0?
    P.Tmin      = 273.15-40;
    P.rho_l     = 1000;
    P.Patmo     = cI.pf - P.rho_l*P.g*cI.Zw; % Pressure at air surface
    P.rho_bg0   = density(cO.pg(end),pI.T0); % pyroclast bubble gas density
    P.LjC       = 1.85; % Jet entrance length scale constant (Turner 1966)
    
    % Setup initial conditions
    n0      = pI.n_0;
    z0      = pI.Ld;
    pa0     = cI.pf - P.rho_l*P.g*pI.Ld;
    W0      = waterProps(pa0,pI.T0);
    rho_s   = (sum(pI.nsi./pI.rhoi))^-1; % Solids bulk density over GSD
    rho_B0  = (pI.n_0/density(pa0,pI.T0) + (1-pI.n_0)./rho_s).^(-1);
%     rho_B0  = (n0/W0.rho + (1-n0)/rho_s)^(-1);
    m_0     = cO.Par.Q;
    m_w0    = m_0*n0;
    m_s0    = m_0*(1-n0);
    psi0    = m_0*pI.u_0;
    numGS   = length(pI.Rgsd);
    
    % Enthalpy
%     C_sb0   = W0.C_p.*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi); % Bulk pyroclast heat capacity
%     H0      = n0*W0.h + (1-n0)*(C_sb0)*(pI.T0 - P.Tr);
    % Alt Enthalpy - more accurate but small diff
    h_s     = enthalpy(P.rho_bg0,pI.T0)*sum(pI.ni.*pI.nsi) ...
                + pI.C_s*sum((1-pI.ni).*pI.nsi)*(pI.T0 - P.Tr);
    H0      = n0*W0.h + (1-n0)*h_s;
    
    gprime  = (pI.rho_B0 - P.rho_l)./max([pI.rho_B0 P.rho_l]) * P.g;
    E0      = m_0 * (H0 + gprime*z0 + 0.5*pI.u_0^2); % Should this be g'* z' to give positive bouyancy in water?
    
%     SSAfromD = @(ri,rhoi,D,A) (3*A*sum(ri.^(2-D))./sum(rhoi.*ri.^(3-D)));
%     SSA0    = SSAfromD(pI.Rgsd,pI.rhoi,pI.D,pI.saScale);
%     flargh change these from pI to wO (must include in ode fun input)

    if isempty(pI.phiFrag_cutoff) % Particle length scale below which quench fragmentation stops
        wO.phi_cutoff = pI.phiFrag_mu; % - pI.phiFrag_sig;
    else
        wO.phi_cutoff = pI.phiFrag_cutoff;
    end
    wO.phi        = -log2(pI.Rgsd*2000);
    wO.phi_pdf    = normpdf(wO.phi,pI.phiFrag_mu,pI.phiFrag_sig);
    wO.phi_pdf    = wO.phi_pdf./sum(wO.phi_pdf); % Ensure sum to 1         % Output PSD from thermal granulation
    wO.fragPari   = pI.fragPar.*(1-wO.phi_pdf./max(wO.phi_pdf));
    wO.fragPari(wO.phi>wO.phi_cutoff) = 0;                                 % Fragmentation energy partition disappears below critical size
    wO.SSAi       = 3*pI.saScale./(pI.rhoi.*pI.Rgsd);                      % Surface area per mass by particle size
    wO.SSA_out    = sum(wO.SSAi.*wO.phi_pdf);                              % Output grain sizes and surface areas are fixed at input for now
    wO.Ess_out    = wO.SSA_out * pI.Es;                                    % Specific surface energy of output grain sizes [J/kg]

    % A0 - Total particle surface area. SSA should be enough?

    % u0 = pI.u_0;
    
    P.Ri0       = P.g*((P.rho_l-pI.rho_B0)/P.rho_l)*pI.r_0/(pI.u_0^2);
    P.De        = 2*pI.r_0;
    P.Lme       = pi^0.25*pI.r_0*abs(P.Ri0)^(-0.5);
    
% ---- Turbulence development scaling ----
    % v1: Jet entrance scaling
%     gprime = P.g*(rho_wa-rho_B)/max([rho_wa rho_B]);
%     Fo =  pi*gprime*(pI.r_0).^2*pI.u_0;
%     Lj =  P.LjC * Mo^(3/4) * Fo^(-1/2);
    
    % v2: Jet length v2 (broken)
%     Mo =  pi*(pI.r_0).^2.*pI.u_0.^2;
%     P.Lj        = Mo^(1/4)*(pI.u_0/pI.r_0).^(1/2);

    % v3: Simple Kotsovinos regime 1
%     P.Lj = (pi/4)^(1/4)*(1/2)^(1/2)*pI.r_0;

    % v4: Bouyancy-Interia balance for start of an MTT plume
    V0 = cO.Par.Q./rho_B0;
    P.Lj = ( pI.u_0 .* pi.^(5/6) .* (pI.r_0.^5./(abs(gprime).*V0)).^(1/3) ).^(3/4);
    
%% SOLVE ODE

    % ============ intitial run ===========
    zSpan = [pI.Ld cI.Zw]; % Z reference frame: Ld = 0? vent = 0?
    IC    = [m_w0 m_s0 psi0 E0 pI.nsi'];
%     nn = zeros(size(IC));
    nn = (5:5+numGS-1);
    options = odeset('RelTol',1e-6,'AbsTol',1e-4,'NonNegative',nn,...
        'Events',@(z,y) stopEntrain1(y,P.xv_crit));
%     tic
    odeSolution = ode15s(@(z,y) waterEntrain(z,y,cI,pI,wO,P),zSpan,IC,options);
%     toc
    
    Z = odeSolution.x';
    Y = odeSolution.y';
    
    % Check stop conditions
%     if (xv < P.xv_crit) && (Z(end) < cI.Zw)
    if (Z(end) < cI.Zw)
        wO.failedPlume = true;
    else
        wO.failedPlume = false;
        
%         % OPTION: Very close to all liquid water fraction, relax tolerances to allow
%         % very rapid density change. Error will increase, but this is
%         % essentially a failed plume anyway.
%         rtol = 1e-3;
%         atol = 1e-3; %[1 10 10 1];
%         refine = 4;
%         options = odeset('RelTol',rtol,'AbsTol',atol,'Refine',refine,'Events',@(z,y) stopEntrain2(y,m_s0));
%         odeSol2 = odextend(odeSolution,@(z,y) waterEntrain(z,y,cI,pI,P),zSpan(2),[],options);
%         % Append solution
    end

    % Y out 
    m_w = Y(:,1);
    wO.m_s = Y(:,2); %m_s0; % COULD include a sedimentation parameterization in future?
    wO.psi = Y(:,3);
    wO.E   = Y(:,4);
    wO.nsi   = Y(:,5:5+numGS-1);

    % ---- Post-process -----
    pI2 = pI;
    
    wO.z = Z;
    wO.m = wO.m_s + m_w;
    wO.u = wO.psi./wO.m;
    wO.H = wO.E./wO.m;
    wO.P = P.rho_l.*P.g.*(cI.Zw - Z) + P.Patmo;
    wO.SSA = sum(wO.nsi.*wO.SSAi',2);
    
    wO.T = zeros(size(Z));
    wO.Tsat  = wO.T;
    wO.rho_w = wO.T;
    wO.xv    = wO.T;
    wO.hw    = wO.T;
    wO.alpha = wO.T;
    wO.r_mix = wO.T;
    for zi = 1:length(Z)
        [~,W] = waterEntrain(Z(zi),Y(zi,:)',cI,pI,wO,P);
        wO.T(zi)     = W.T;
        wO.Tsat(zi)  = W.Tsat;
        wO.rho_w(zi) = W.rho_w;
        wO.xv(zi)    = W.xv;
        wO.hw(zi)    = W.hw;
        wO.alpha(zi) = W.alpha;
        wO.r_mix(zi) = W.r_mix;
    end
%     toc
    wO.Lj   = P.Lj;
    wO.m_v = m_w.*wO.xv;
    wO.m_l = m_w.*(1-wO.xv);
    wO.n_v = wO.m_v./wO.m;
    wO.n_l = wO.m_l./wO.m;
    n_w = m_w./wO.m;
    wO.rho_B = (n_w./wO.rho_w + (1-n_w)./rho_s).^(-1);
    wO.r = (wO.m./(pi.*wO.rho_B.*wO.u)).^(1/2);

    % Generate new PSD - single step at the end for now
%     D_all   = (2.8:.02:4);
%     SSA_all = SSAfromD(pI.Rgsd,pI.rhoi,D_all,pI.saScale);
%     wO.D    = interp1(SSA_all,D_all,wO.A(end),'spline');
%     pI2.D   = wO.D;
%     [pI2.nsi,pI2.ni,pI2.Rgsd,pI2.rhoi,pI2.porR] = getPSD(cI.rho_melt,pI.phi0, cO.pg(end), wO.T(end), wO.xv(end), wO.D,pI.phiSz_min,pI.phiSz_max);

    
    % Update plume params    
    pI2.n_0     = (m_w(end) + wO.m_s(end).*(sum(pI.ni.*pI.nsi) - sum(pI2.ni.*pI2.nsi))) ./ (m_w(end) + wO.m_s(end));
    pI2.T0      = wO.T(end);
    pI2.vh0     = cI.vh0 + wO.z(end);
    pI2.u_0     = wO.u(end);
    pI2.rho_B0  = wO.rho_B(end);
    pI2.r_0     = ((wO.m_s(end) + m_w(end))./(pi*wO.rho_B(end)*wO.u(end))).^(1/2); % PROBLEM WITH PLUME RADIAL EXPANSION THAT IS FAR TOO RAPID HERE
    pI2.xv_0    = wO.xv(end);
    pI2.nsi     = wO.nsi(end,:)';
%     pI2.rho_g0  = 
    % m_w(end) + dm_ws;
    
    % Run enthalpy calcs for m_v, m_l?
% pI.T0
%     toc

end

%% Functions

function dT = findT(T,P,H,n,n_gs,C_s,rho_gb,Tr)
% Iteration function to find T from total enthalpy
%   T       = guess temp
%   P       = pressure
%   H       = actual enthalpy
%   n       = free gas mass fraction
%   n_gs    = bubble gas mass fraction (relative to total pyroclast mass)
%   C_s     = melt heat capacity
%   rho_gb  = bubble gas density
%   Tr      = reference Temperature

%     W = waterProps(P,T);
%     rho_w = density(P,T);    % kg/m^3
    hw    = enthalpy(density(P,T),T); % J/kg
    
    [h_bg,~] = enthalpy(rho_gb,T);
    h_s = h_bg*n_gs + C_s*(1-n_gs)*(T - Tr); 
    
    H_est = (1-n).*h_s + n.*hw;
    dT = abs(H - H_est);
end

% function dE = findT_2(T,P,Tm,hw0,C_s,F,dm_w,m_s)
%     
%     hw = enthalpy(density(P,T),T);
%     dE = abs(dm_w*(hw - hw0) + m_s*F*C_s*(T-Tm));
% %     left = F*C_s*Tm - dm_w/m_s*hw0;
% %     right = F*C_s
% end

function [dy_dz WO] = waterEntrain(z,y,cI,pI,wO,P)
% OPTIONS: Gas release at each frag step? Requires GSD calc...
%
%

persistent rho_prev %ax % zchk % hw_prev Tprev dm_w_pr m_w_pr ax z_pr
global xv
if isempty(rho_prev)
    rho_prev = pI.rho_B0;
%     zchk = false;
%     tic
%     figure(2)
%     ax(1)=subplot(2,1,1);
%     ax(2)=subplot(2,1,2);
%     ax(3)=subplot(3,1,3);

end
%     zvals = [70:10:170];
%     zchk = false;
%     if any(abs(zvals-z)<1) && ~zchk
%         tic
%     end

    numGS = length(pI.Rgsd);
%     z = y(1);
    m_w = y(1);
    m_s = y(2); %cI.Q*(1-pI.n_0);
    psi = y(3);
    E   = y(4);
    nsi = y(5:5+numGS-1);
%     A   = y(4);
    
    m = m_s + m_w;
    n = m_w/m;
    u = psi/m;
    pa = P.rho_l.*P.g.*(cI.Zw - z) + P.Patmo;  % Ambient pressure
    rho_wa  = density(pa,pI.Tw0);
    gprime = (rho_prev - rho_wa)./max([rho_prev rho_wa]) * P.g;
    H = E/m - u^2/2 - gprime*z;

    % Particle density/gas frac assuming bubble gas at equal pressure to ambient
%     rho_w_prev = density(pa,pI.T0);
%     C      = pI.pori./(1-pI.pori).*rho_w_prev/pI.rho_m;
%     ni_p   = C./(1+C);
%     rho_i  = (1-pI.pori)*pI.rho_m + (pI.pori)*rho_w_prev; 
%     rho_s  = 1./(sum(pI.nsi./rho_i));
     % Density change is very small above Tsat...

    % ----- Get water properties and mixture Temperature ------
    W0 = waterProps(pa,pI.T0); 
%     C_sb = W0.C_p.*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi); % Bulk pyroclast heat capacity
    Tmax = H/((1-n)*pI.C_s) + P.Tr; % Use C_b here?
%     Tmax = min([H/((1-n)*pI.C_s) + P.Tr pI.T0]); % Use C_b here?

    % Pyroclast enthalpy
    % --> Bubble gas is minor component, neglect it's thermodynamics and 
    %     only conserve mass for now. 
    %     Room here for interesting work on pyroclast water ingestion and 
    %     related effects
%     pb = pressure(P.rho_bg0,W0.Tsat);
    [h_bgS,~] = enthalpy(P.rho_bg0,W0.Tsat);
    h_s_sat = h_bgS*sum(pI.ni.*nsi) ...
            + pI.C_s*sum((1-pI.ni).*nsi)*(W0.Tsat - P.Tr); 
        
    Hvap = (1-n)*h_s_sat + n*W0.h_vs;
    Hliq = (1-n)*h_s_sat + n*W0.h_ws;
    
    Tfun = @(x) findT(x,pa,H,n,sum(pI.ni.*pI.nsi),pI.C_s,P.rho_bg0,P.Tr);
    opts.TolX = 1e-3;
    
    if H > Hvap                     % Water is vapour
        T = fminbnd(Tfun,W0.Tsat,Tmax,opts);
        rho_w = density(pa,T);
        xv = 1;
    elseif and(H<=Hvap, H>=Hliq)    % Water is saturated
        xv = (H - Hliq) / (Hvap-Hliq); % Dryness fraction
        T = W0.Tsat;
        rho_w = densityTX(W0.Tsat,xv);
    else                            % Water is liquid
        T = fminbnd(Tfun,P.Tmin,W0.Tsat,opts);
        rho_w = density(pa,T);
        xv = 0;
    end
    
    % ---------------------------------------------------------

    
    hw      = enthalpy(rho_w,T);                    % Plume water enthalpy
    hw0     = enthalpy(rho_wa,pI.Tw0);              % Ambient water enthalpy
    
    rho_s   = (sum(pI.nsi./pI.rhoi))^-1;            % Solids bulk density (assuming no change in D)
    rho_B   = (n/rho_w + (1-n)/rho_s)^-1;           % Plume bulk density
    r       = (m/(pi*rho_B*u))^(1/2);                  % Plume radius
    
    %% ---------------- Entrainment model ----------------
    Ri      = P.g*((rho_wa-rho_B)/rho_wa)*r/(u^2); % Local Richardson #
    alphaf  = alpha_cara(Ri,z,P.De,P.Lme,r);
    if alphaf<0.05
        alphaf=0.05;
    elseif alphaf>0.17
        alphaf=0.17;
    end
    
    % Rayleigh-Taylor entrainment from Zhang et al 2020
    if pI.alphaRTscale~=0
       c_rt = alphaRT(pI.r_0,pI.Tw0,pa,cI.Q,rho_B,u,r);
       alphaf = (1-pI.alphaRTscale)*alphaf + pI.alphaRTscale*(c_rt);
    end
    % Kelvin-Helmholtz entrainment from Zhang et al 2020
%     if pI.alphaRTscale~=0
%        c_kh = alphaKH(pI.r_0,cI.Q,rho_w,u,r);
%        alphaf = (1-pI.alphaRTscale)*alphaf + pI.alphaRTscale*(c_kh);
%     end    
    
    % Jet Entrance Option 1: Simple linear scale for mixing layer development over a jet entrance length
%     if and( pI.useJetEntranceLength , (z-pI.Ld)<P.Lj)
%         r_mix  = r*(z-pI.Ld)./P.Lj;
%         alphaf = alphaf * (z-pI.Ld)./P.Lj; 
%     else
%         r_mix = r;
%     end
    
    % Jet Entrance Option 2: Non-linear (volumetric) mixing layer/entrainment scale
%     gprime = (rho_B - rho_wa)./max([rho_B rho_wa]) * P.g;
%     if pI.useJetEntranceLength
%         r_mix  = 2*((alphaf^2 * abs((r-pI.r_0)/gprime))^(1/2) * u); % Width of mixing layer
%         if r_mix<0; r_mix = 0; elseif r_mix>r; r_mix = r; end
% %         r_mix_frac = r_mix./r;
%         m_s_frac   = (2*r*r_mix - r_mix^2)./r^2;
%         alphaf     = alphaf * m_s_frac; %r_mix_frac;
% %     else
% %         r_mix_frac = 1;
% %         m_s_frac = 1;
%     end

    % Jet Entrance Option 3: Hybrid - linear jet entrance with volumetric scale
    if and( pI.useJetEntranceLength , (z-pI.Ld)<P.Lj)
        r_mix  = r*(z-pI.Ld)./P.Lj;
        m_s_frac   = (2*r*r_mix - r_mix^2)./r^2;
        alphaf = alphaf * m_s_frac; 
    else
        r_mix = r;
    end   
    
%     figure(2)
%     scatter(ax(1),z,r_mix/r,'ob')
%     hold(ax(1),'on')
%     scatter(ax(2),z,alphaf,'ob')
%     hold(ax(2),'on')
%     pause(0.1)
    
    u_eps = alphaf * u * sqrt(rho_B * rho_wa); % Wrap up the density in u_eps for now
    % ------------------------------
    
    dm_w_ent = 2 * pi * r * u_eps;  % rho = sqrt(rho*rho_a)? % Entrained water
    
    
%%  % -------- Fragmentation and Grain size output model  -------

    % Entrained water enthalpy change estimate
    [dhw,~] = getDeltaHw(hw,hw0,W0.h_ws,W0.h_vs,pI.C_s,m_s,m_w,dm_w_ent,rho_w,T,W0.Tsat,pI.fragPar);
    % --------
    
    % Glass transition limit on fragmentation energy
    hs_sm = (@(x) 1./(1 + exp(-2*(3./pI.T_g_rng)*(x-(pI.T_g+pI.T_g_rng/2)))));
    fragPari = wO.fragPari*hs_sm(T);
%     fragPari(nsi<0.0001)=0;
    
%     if T > (pI.T_g + pI.T_g_rng)
%         fragPari = pI.fragPari;
%     elseif and(T >= pI.T_g, T <= (pI.T_g + pI.T_g_rng))
%         fragPari = pI.fragPari * (T - pI.T_g)/(pI.T_g_rng);
%     else
%         fragPari = 0*pI.fragPar;
%     end
%     fragPari    = pI.fragPari; % or not

%     fragPar_eff = sum((1-pI.ni).*pI.nsi.*fragPari);          % Effective fragmentation partitioning, accounting for solids PSD
    fragPar_eff = sum((1-pI.ni).*nsi.*fragPari);          % Effective fragmentation partitioning, accounting for solids PSD
%     fragPar_eff(nsi<0.0001) = 0;
    
     % Compare m_frag (v3) w/ dE_ss below
    dE_ss = fragPar_eff/(1-fragPar_eff)*(dhw) * dm_w_ent/(m_s); % Specific fragmentation energy, J/kg_melt
    
    % Distribution of mass removed from coarse fraction
    nsi_in  = nsi.*fragPari.*wO.SSAi.*(1-pI.ni)./sum(nsi.*fragPari.*wO.SSAi.*(1-pI.ni)); % Scaled by available surface area
%     nsi_in  = nsi.*fragPari.*(1-pI.ni)./sum(nsi.*fragPari.*(1-pI.ni));                   % Scaled by available mass fraction
    
    dn_s_frag = dE_ss/wO.Ess_out; % Total solids mass fraction fragmented at this step

    %% Deltas
    gprime = (rho_B - rho_wa)./max([rho_B P.rho_l]) * P.g; % Reduced grav for bouyancy force and potential E
    
    dnsi_dz = dn_s_frag * (-nsi_in + wO.phi_pdf);            % PSD change
    dm_w_frag = m_s*dn_s_frag*(sum((pI.ni).*(nsi)./(1-pI.ni)) - sum((pI.ni).*(nsi+dnsi_dz)./(1-pI.ni))); % Gas release from fragmentation

    dm_w_dz = dm_w_ent + dm_w_frag;                         % Water mass
    dm_s_dz = -dm_w_frag;                                   % Solids gas release
    dpsi_dz = P.g*(rho_wa - rho_B)*r^2;                    % Momentum
    dE_dz   = dm_w_ent * (gprime*z + hw0) - m_s*dE_ss;        % Total enthalpy change
%     dA_dz = dE_ss / pI.Es;                                % Specific surface area, m^2/kg_melt/m

    
    % Incremental fragmentation and grain size change COULD happen here
    %   --> Allows incremental gas release from pyroclasts (generally
    %       likely to be small at this stage, but could be significant
    %       if introduced all at once later)
    %   --> Allows a parameterization for solids mass loss later
    %   --> Requires initial calculation of entrained water, then second
    %   calculation of d_total_water = w_entrain + w_frag_release
    
    dy_dz    = zeros(4+numGS,1);
    dy_dz(1) = dm_w_dz;
    dy_dz(2) = dm_s_dz;
    dy_dz(3) = dpsi_dz;
    dy_dz(4) = dE_dz;
    dy_dz(5:5+numGS-1) = dnsi_dz;
%     dy_dz(5) = dA_dz;
    
    rho_prev = rho_B;
%     hw_prev  = hw;
%     Tprev    = T;
%     dm_w_pr  = dm_w_dz;
%     m_w_pr   = m_w;
%     z_pr     = z;
    if nargout>1
        WO.T     = T;
        WO.Tsat  = W0.Tsat;
        WO.rho_w = rho_w;
        WO.xv     = xv;
        WO.hw    = hw;
        WO.alpha = alphaf;
        WO.r_mix = r_mix;
    end
    
    if or( any(isnan(y)) , isnan(rho_prev) )
        fprintf('Z: %.2f, H: %.2e, Tsat: %.2f, T: %.2f, rho_prev: %.2f\n',z,H,W0.Tsat,T,rho_prev)
    end
%     if any(abs(zvals-z)<1) && ~zchk
%         fprintf('%.2f, %.2f\n',z,toc)
%         zchk = true;
%         tic
%     end
%     if ~any(abs(zvals-z)<1) && zchk
%         zchk = false;
%     end

end

function [value,isterminal,direction] = stopEntrain1(y,xv_crit)
    % Stop integration of velocity reaches zero, or water dryness fraction
    % gets close to 0

    global xv
    
    % u = 0, T ~ Tw0?
    m_w = y(1);
    m_s = y(2);
    psi = y(3);
    m   = m_w + m_s;
    u   = psi/m;

    value      = [u xv-xv_crit];
    isterminal = [1 1];
    direction  = [0 -1];

end

% function [value,isterminal,direction] = stopEntrain2(y,m_s)
%     % Stop integration if velocity reaches zero only
% 
%     
%     % u = 0, T ~ Tw0?
%     m_w = y(2);
%     psi = y(3);
%     m   = m_w + m_s;
%     u   = psi/m;
% 
%     value      = u;
%     isterminal = 1;
%     direction  = 0;
% 
% end

function [dhw,Tf,x] = getDeltaHw(hw,hw0,h_ws,h_vs,c_s,m_s,m_w,dm_w,rho_w,T,Tsat,f)
% Approximation for change in water enthalpy as function of water
% entrained. 
%   --> For calculation of fragmentation energy only.
%   --> Only handles temperature decrease, not increase.
%   dhw = change in enthalpy of ambient water on entrainment
%   Tf  = approximate final temperature of mixture at equilibrium
%   x   = dryness fraction of final mixture

    if hw > h_vs
        cpv = heatCapacityP(rho_w,T);
%         Tf  = -(hw - hw0) / ( (1-f)*m_s/dm_w*c_s - cpv) + T;
%         gprime = (rho_w - 1000)./max([rho_w 1000]) * 9.81;
        Tf = ((1-f)*m_s*c_s*T + m_w*cpv*T - dm_w*(hw - cpv*T - hw0)) / (cpv*(m_w+dm_w) + (1-f)*m_s*c_s);
        
        dhw = hw + cpv*(Tf - T) - hw0;
        
        if hw0+dhw < h_vs
            if T<Tsat; error('Somehow T got real low...'); end
            dh1 = ((1-f)*c_s*(T-Tsat)*m_s + m_w*(hw-h_vs)) / dm_w;
            [dhw2,Tf,x] = getDeltaHw(h_vs,hw0+dh1,h_ws,h_vs,c_s,m_s,m_w,dm_w,rho_w,Tsat,Tsat,f);
            dhw      = dhw2 + dh1;
%             x   = (hw0+dhw - h_ws) / ();
        else
            x = 1;
        end
        
    elseif and( h_ws<= hw, hw <= h_vs)
        dhw = (m_w*hw + dm_w*hw0) / (m_w + dm_w) - hw0;
        Tf = Tsat;
        x   = (hw0 + dhw - h_ws) / (h_vs - h_ws);
        
        if hw0+dhw < h_ws
            dh1 = (hw - h_ws) * m_w / dm_w;
            [dhw2,Tf,x] =  getDeltaHw(h_ws-1e-3,hw0+dh1,h_ws,h_vs,c_s,m_s,m_w,dm_w,rho_w,Tsat,Tsat,f);
            dhw      = dhw2 + dh1;
%             dhw      = (hw - h_ws) + dhw;
        end
        
    else
        cpw = heatCapacityP(rho_w,T);
%         Tf  = -(hw - hw0) / ( (1-f)*m_s/dm_w*c_s - cpw) + T;
        Tf = ((1-f)*m_s*c_s*T + m_w*cpw*T - dm_w*(hw - cpw*T - hw0)) / (cpw*(m_w+dm_w) + (1-f)*m_s*c_s);
        dhw = hw + cpw*(Tf - T) - hw0;
        x = 0;
    end
    
end

% function SSA = SSAfromD(ri,rhoi,D,A)
% % Calculate specific surface area, m^2/kg, from power-law grain size distribution
% %   ri   = vector of effective particle radii
% %   rhoi = density of particles at each radii
% %   D    = power law exponent
% %   A    = roughness scale parameter
% 
%     SSA = (3*A*sum(ri.^(2-D))./sum(rhoi.*ri.^(3-D)));
% 
% end