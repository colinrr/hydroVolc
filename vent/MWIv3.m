function [pI2,wO] = MWIv3(cI,cO,pI,dO)
% IN: 
%       pI    plume input params
%       dO    decompression ouput

%     P.rho_m   = cI.rho_melt;      % density of solid pyroclasts [kg/m^3]
%     P.n_0     = pI.n_0;           % exit gas content            [wt.%]
%     P.r_0     = pI.r_0;           % jet radius                  [m]
%     P.u_0     = pI.u_0;           % exit velocity               [m/s]
%     P.T0      = pI.T0;            % Magma initial temperature   [K]
%     P.Tw0     = 273.15 + 1;       % Water initial temperature   [k]
%     P.g       = 9.81;             % gravity (m/s^2)
%     Tref
% 
%     Es      = 100;                % surface energy of pyroclasts  [J/m^2]
%     Cpm     = 1250;               % magma spec. heat capacity     [J/(kgK)]
    
    % Need:
    %  10  Roughness scale
    %  0.1 Fragmentation energy partitioning fraction
    
%     Cpv     = 1952;             % spec. heat of volcanic gas (water) at constant pressure (J/kg/K)
%     Cpl     = 4190;             % spec. heat of liquid water      (J/kg/K)
%     R_v     = 461;              % gas constant of volcanic gas (water) [J/kg/K]
%     L       = 2.257e6;          % latent heat of vaporization (J kg^-1)
%     rho_l   = 1000;             % density liquid water (kg/m^3)
    
%     dTmin = 150;   % Minimum temperature drop (very ish)

%     fragPar = 1;

%% SETUP
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
%     rho_B0  = (n0/W0.rho + (1-n0)/rho_s)^(-1);
    m_0     = cI.Q;
    m_w0    = m_0*n0;
    m_s0    = m_0*(1-n0);
    psi0    = m_0*pI.u_0;
    numGS   = length(pI.Rgsd);
    
    % Enthalpy
%     C_sb0   = W0.C_p.*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi); % Bulk pyroclast heat capacity
%     H0      = n0*W0.h + (1-n0)*(C_sb0)*(pI.T0 - P.Tr);
    % Alt Enthalpy - more accurate but small diff
%     h_s     = enthalpy(P.rho_bg0,pI.T0)*sum(pI.ni.*pI.nsi) ...
%                 + pI.C_s*sum((1-pI.ni).*pI.nsi)*(pI.T0 - P.Tr);
    h_s     = pI.C_s*(pI.T0 - P.Tr);
            
    % Original energy calc
%     H0      = n0*W0.h + (1-n0)*h_s;
%     gprime  = (pI.rho_B0 - P.rho_l)./max([pI.rho_B0 P.rho_l]) * P.g;
%     E0      = m_0 * (H0 + gprime*z0 + 0.5*pI.u_0^2); % Should this be g'* z' to give positive bouyancy in water?
    
    % New energy calc (thermal disequilibrium between particles and water)
    gprime_w = (W0.rho - P.rho_l)./max([W0.rho  P.rho_l]) * P.g;
    gprime_s = (rho_s - P.rho_l)./max([W0.rho  P.rho_l]) * P.g;
    Ew0     = m_w0 * (W0.h + gprime_w*z0 + 0.5*pI.u_0^2);
    Es0     = m_s0  * (h_s  + gprime_s*z0 + 0.5*pI.u_0^2); 
    
%     SSAfromD = @(ri,rhoi,D,A) (3*A*sum(ri.^(2-D))./sum(rhoi.*ri.^(3-D))); % Relies on power-law psd
    SSAfromD = @(ri,rhoi,nsi) sum((3*pI.saScale./(rhoi.*ri)).*nsi);
    SSA0    = SSAfromD(pI.Rgsd,pI.rhoi,pI.nsi);
    % A0 - Total particle surface area. SSA should be enough?
    
    % Params for output particle distribution from thermal granulation
    flargh change these from pI to wO
    pI.phi_cutoff = pI.phiFrag_mu - pI.phiFrag_sig;                        % Particle length scale below which quench fragmentation stops
    pI.phi        = -log2(pI.Rgsd*2000);
    pI.phi_pdf    = normpdf(pI.phi,pI.phiFrag_mu,pI.phiFrag_sig);
    pI.phi_pdf    = pI.phi_pdf./sum(pI.phi_pdf); % Ensure sum to 1         % Output PSD from thermal granulation
    pI.fragPari   = pI.fragPar.*(1-pI.phi_pdf./max(pI.phi_pdf));
    pI.fragPari(pI.phi>pI.phi_cutoff) = 0;                                 % Fragmentation energy partition disappears below critical size
    pI.r_crit     = 2.^(-pI.phiFrag_mu)/1000;                              % Conductivity length scale corresponds to mean output grain size
    pI.SSAi       = 3*pI.saScale./(pI.rhoi.*pI.Rgsd);                      % Surface area per mass by particle size
    pI.SSA_out    = sum(pI.SSAi.*pI.phi_pdf);                              % Output grain sizes and surface areas are fixed at input for now
    pI.Ess_out    = pI.SSA_out * pI.Es;                                    % Specific surface energy of output grain sizes [J/kg]

    % u0 = pI.u_0;
    
    P.Ri0       = P.g*((P.rho_l-pI.rho_B0)/P.rho_l)*pI.r_0/(pI.u_0^2);
    P.De        = 2*pI.r_0;
    P.Lme       = pi^0.25*pI.r_0*abs(P.Ri0)^(-0.5);
    
    % Jet entrance scaling
%     gprime = P.g*(rho_wa-rho_B)/max([rho_wa rho_B]);
%     Mo =  pi*(pI.r_0).^2.*pI.u_0.^2;
%     Fo =  pi*gprime*(pI.r_0).^2*pI.u_0;
%     Lj =  P.LjC * Mo^(3/4) * Fo^(-1/2);
%     P.Lj        = Mo^(1/4)*(pI.u_0/pI.r_0).^(1/2);
    
%% SOLVE ODE
    abstol = [1e-3]; % 1e-6
    options = odeset('RelTol',1e-5,'AbsTol',abstol,'Events',@(z,y) stopEntrain(y,m_s0));
    zSpan = [pI.Ld cI.Zw]; % Z reference frame: Ld = 0? vent = 0?
    % Initial conditions: m_w0 Em0 Ew0 psi0 A nsi?
    IC    = [m_w0 m_s0 psi0 Ew0 Es0 pI.nsi'];
    
    tic
    odeSolution = ode15s(@(z,y) waterEntrain(z,y,cI,pI,P),zSpan,IC,options);
    toc
    Z = odeSolution.x';
    Y = odeSolution.y';

    % Y out 
    m_w    = Y(:,1);
    m_s    = m_s0; % COULD include a sedimentation parameterization in future?
    wO.psi = Y(:,2);
    wO.Ew  = Y(:,3);
    w0.Em  = Y(:,4);
    wO.A   = Y(:,5);

    % ---- Post-process -----
    pI2 = pI;
    
    wO.z = Z;
    wO.m = m_s + m_w;
    wO.u = wO.psi./wO.m;
    wO.Hw = wO.Ew./wO.m_w;
    wO.Hs = m_s./wO.Em;
    wO.P = P.rho_l.*P.g.*(cI.Zw - Z) + P.Patmo;
    
    % Get T
%     tic
    flargh - check_calcs_below;
    wO.T = zeros(size(Z));
    wO.Tsat = wO.T;
    wO.rho_w = wO.T;
    wO.xv     = wO.T;
    wO.hw    = wO.T;
    wO.alpha = wO.T;
    for zi = 1:length(Z)
        [~,W] = waterEntrain(Z(zi),Y(zi,:),cI,pI,P);
        wO.T(zi)     = W.T;
        wO.Tsat(zi)  = W.Tsat;
        wO.rho_w(zi) = W.rho_w;
        wO.xv(zi)     = W.xv;
        wO.hw(zi)    = W.hw;        
        wO.alpha(zi) = W.alpha;
    end
%     toc
    wO.m_v = m_w.*wO.xv;
    wO.m_l = m_w.*(1-wO.xv);
    wO.n_v = wO.m_v./wO.m;
    wO.n_l = wO.m_l./wO.m;
    n_w = m_w./wO.m;
    wO.rho_B = (n_w./wO.rho_w + (1-n_w)./rho_s).^(-1);
    wO.r = (wO.m./(pi.*wO.rho_B.*wO.u)).^(1/2);

    % Generate new PSD - single step at the end for now
    D_all   = (2.8:.02:4);
    SSA_all = SSAfromD(pI.Rgsd,pI.rhoi,D_all,pI.saScale);
    wO.D    = interp1(SSA_all,D_all,wO.A(end),'spline');
    pI2.D   = wO.D;
    [pI2.nsi,pI2.ni,pI2.Rgsd,pI2.rhoi,pI2.pori] = getPSD(cI.rho_melt,pI.phi0, cO.pg(end), wO.T(end), wO.D,pI.phiSz_min,pI.phiSz_max);

    % Update plume params    
    pI2.n_0     = (m_w(end) + m_s.*(sum(pI.ni.*pI.nsi) - sum(pI2.ni.*pI2.nsi))) ./ (m_w(end) + m_s);
    pI2.T0      = wO.T(end);
    pI2.vh0     = cI.vh0 + wO.z(end);
    pI2.u_0     = wO.u(end);
    pI2.rho_B0  = wO.rho_B(1);
    pI2.r_0     = ((m_s + m_w(end))./(pi*wO.rho_B(end)*pI.u_0)).^(1/2); % PROBLEM WITH PLUME RADIAL EXPANSION THAT IS FAR TOO RAPID HERE
    
    % m_w(end) + dm_ws;
    
    % Run enthalpy calcs for m_v, m_l?
% pI.T0
    toc

end

%% Functions

function dT = findT_old(T,P,H,n,n_gs,C_s,rho_gb,Tr)
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
    h_s_sat = h_bg*n_gs + C_s*(1-n_gs)*(T - Tr); 
    
    H_est = (1-n).*h_s_sat + n.*hw;
    dT = abs(H - H_est);
end

function dT = findTpyroclast(T,P,H,n_gs,C_s,rho_gb,Tr)
% Iteration function to find T of pyroclasts from their total enthalpy,
% assuming bubble gas is in pressure equilibrium with exterior and
% temperature equilibrium with pyroclast.
%   T       = guess temp
%   P       = pressure
%   H       = actual enthalpy
%   Tr      = reference Temperature

%     W = waterProps(P,T);
%     rho_w = density(P,T);    % kg/m^3
    h_bg    = enthalpy(density(P,T),T); % J/kg
    
%     [h_bg,~] = enthalpy(rho_gb,T);
    h_s_sat = h_bg*n_gs + C_s*(1-n_gs)*(T - Tr); 
    
%     H_est = (1-n).*h_s_sat + n.*hw;
    dT = abs(H - h_s_sat);
end

function dT = findTwater(T,P,H)
% Iteration function to find T of pyroclasts from their total enthalpy,
% assuming bubble gas is in pressure equilibrium with exterior and
% temperature equilibrium with pyroclast.
%   T       = guess temp
%   P       = pressure
%   H       = actual enthalpy
%   n_gs    = bubble gas mass fraction (relative to total pyroclast mass)
%   C_s     = melt heat capacity
%   rho_gb  = bubble gas density
%   Tr      = reference Temperature

%     W = waterProps(P,T);
%     rho_w = density(P,T);    % kg/m^3
    h_w    = enthalpy(density(P,T),T); % J/kg
    
%     [h_bg,~] = enthalpy(rho_gb,T);
%     h_s_sat = h_bg*n_gs + C_s*(1-n_gs)*(T - Tr); 
    
%     H_est = (1-n).*h_s_sat + n.*hw;
    dT = abs(H - h_w);
end

function dT = findTequil(T,P,m_w,hwi,m_s,C_s,T_s,fragPar)
% Find enthalpy change that would occur for water to reach equilibrium with
% particles
hwf = enthalpy(density(P,T),T);
T_est = m_w * (hwf - hwi) / ((1-fragPar)*C_s.*m_s) + T_s;

dT = abs(T_est - T);

end
% function dE = findT_2(T,P,Tm,hw0,C_s,F,dm_w,m_s)
%     
%     hw = enthalpy(density(P,T),T);
%     dE = abs(dm_w*(hw - hw0) + m_s*F*C_s*(T-Tm));
% %     left = F*C_s*Tm - dm_w/m_s*hw0;
% %     right = F*C_s
% end

function [dy_dz WO] = waterEntrain(z,y,cI,pI,P)
% OPTIONS: Gas release at each frag step? Requires GSD calc...
%
%
%     tic
    persistent rho_w_prev ax % hw_prev Tprev dm_w_pr m_w_pr ax z_pr
    if isempty(rho_w_prev)
        rho_w_prev = pI.rho_g0;
        
%         ax(1) = subplot(3,1,1);
%         ax(2) = subplot(3,1,2);
%         ax(3) = subplot(3,1,3);
%         hold(ax,'on')
    end
    numGS = length(pI.Rgsd);
    
%     z = y(1);
%     m_w = y(1);
%     m_s = cI.Q*(1-pI.n_0);
%     psi = y(2);
%     Ew  = y(3);
%     Es  = y(4);
%     A   = y(5);

    % dm_w_dz, dm_s_dz, dQ_s_dz, dQ_w_dz, dnsi_dz, dpsi_dz
    
    m_w = y(1);
    m_s = y(2);
    psi = y(3);
    Ew  = y(4);
    Es  = y(5);
    nsi = y(6:6+numGS-1);
    
    
    m = m_s + m_w;
    n = m_w/m;
    u = psi/m;
    pa = P.rho_l.*P.g.*(cI.Zw - z) + P.Patmo;  % Ambient pressure
    
    % Old E calc
%     gprime = (rho_prev - P.rho_l)./max([rho_prev P.rho_l]) * P.g;
%     H = E/m - u^2/2 - gprime*z;
    
    % Particle density/gas frac assuming bubble gas at equal pressure to ambient
    C      = pI.pori./(1-pI.pori).*rho_w_prev/pI.rho_m;
    ni_p   = C./(1+C);
    rho_i  = (1-pI.pori)*pI.rho_m + (pI.pori)*rho_w_prev; 
    rho_s  = 1./(sum(nsi./rho_i));
    
    rho_wa  = density(pa,pI.Tw0); % Ambient water
    % New disequilibrium E calc    
    gprime_w = (rho_w_prev - rho_wa)./max([rho_w_prev  P.rho_l]) * P.g; % Need persistent
    gprime_s = (rho_s - rho_wa)./max([rho_s  P.rho_l]) * P.g; % Need persistent
%     gprime = (rho_prev - P.rho_l)./max([rho_prev P.rho_l]) * P.g;
    Hs = Es/m_s - u^2/2 - gprime_s*z;
    Hw = Ew/m_w - u^2/2 - gprime_w*z;

    % ----- Get water properties and mixture Temperature ------
%     Tw = temperature?         % Consider using persistent T_prev here?
    W0 = waterProps(pa,pI.T0);  % Consider using persistent T_prev here?
    
    
%     C_sb = W0.C_p.*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi); % Bulk pyroclast heat capacity
%     Tmax_s = Hs/((1-n)*pI.C_s) + P.Tr; % Use C_b here?
    Tmax_w = temperature(pa,W0.rho);
%     Tmax = min([H/((1-n)*pI.C_s) + P.Tr pI.T0]); % Use C_b here?

%     C_sb   = W0.C_p.*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi)
%     Hvap = (1-n)*C_sb*(W0.Tsat - P.Tr) + n*W0.h_vs;
%     Hliq = (1-n)*C_sb*(W0.Tsat - P.Tr) + n*W0.h_ws;

    % Pyroclast enthalpy
    % --> Bubble gas is minor component at this point, assume pressure balanced with 
    %     ambient for now and neglect it's enthalpy (order <1%)
    %     Room here for interesting work on pyroclast water ingestion and related effects
%     pb = pressure(P.rho_bg0,W0.Tsat);

%     [h_bgS,~] = enthalpy(W0.rho_vs,W0.Tsat);
%     h_s_sat = pI.C_s*sum((1-pI.ni).*pI.nsi)*(W0.Tsat - P.Tr);
    T_s = Hs./(pI.C_s) + P.Tr; % Solids temp
    try
        rho_bw = density(pa,T_s);
    catch
        disp('Whelp, shit')
    end
    H_bw = enthalpy(rho_bw,T_s); % Bubble gas
        
%     Hvap = (1-n)*h_s_sat + n*W0.h_vs;
%     Hliq = (1-n)*h_s_sat + n*W0.h_ws;
%     toc
    
    Tfun = @(x) findTwater(x,pa,Hw);
    opts.TolX = 1e-3;
    if Hw > W0.h_vs                     % Water is vapour
        T_w = fminbnd(Tfun,W0.Tsat,Tmax_w,opts);
        rho_w = density(pa,T_w);
        xv = 1;
    elseif and(Hw<=W0.h_vs, Hw>=W0.h_ws)    % Water is saturated
        xv = (Hw - W0.h_ws) / (W0.h_vs-W0.h_ws);
        T_w = W0.Tsat;
        rho_w = densityTX(W0.Tsat,xv);
    else                            % Water is liquid
        T_w = fminbnd(Tfun,P.Tmin,W0.Tsat,opts);
        rho_w = density(pa,T_w);
        xv = 0;
    end
      % !!!NOTE!!!: Possible that the linear transition to xv=0 near
      % h_ws is causing the integrator to fail. Possibly try a gaussian
      % smooth here?
    % ---------------------------------------------------------

%     toc
%     hw      = enthalpy(rho_w,T);                    % Plume water enthalpy
    hw0     = enthalpy(rho_wa,pI.Tw0);              % Ambient water enthalpy
    
%     rho_s   = (sum(pI.nsi./pI.rhoi))^-1;            % Solids bulk density (assuming no change in D)
    rho_B   = (n/rho_w + (1-n)/rho_s)^-1;           % Plume bulk density
    r       = (m/(pi*rho_B*u))^(1/2);                  % Plume radius
    
    % Equilibrate bubbles
    C      = pI.pori./(1-pI.pori).*rho_w/pI.rho_m;
    ni    = C./(1+C);
    dm_w_bub = m_s.*(sum(ni_p.*nsi) - sum(ni.*nsi));
    
    % ---- Entrainment model ------
    Ri      = P.g*((P.rho_l-rho_B)/P.rho_l)*r/(u^2); % Local Richardson #
    alphaf  = alpha_cara(Ri,z,P.De,P.Lme,r);
    if alphaf<0.05
        alphaf=0.05;
    elseif alphaf>0.17
        alphaf=0.17;
    end
    
    % Simple scale for mixing layer development over a jet entrance length
%     if and( pI.useJetEntranceLength , (z-pI.Ld)<P.Lj)
    gprime = (rho_B - rho_wa)./max([rho_B rho_wa]) * P.g;
    if pI.useJetEntranceLength    
        r_mix = min([r  2*((alphaf^2 * abs((r-pI.r_0)/gprime))^(1/2) * u)]);  % Width of mixing layer
        if r_mix<0; r_mix = 0; end
        r_mix_frac = r_mix./r;
        m_s_frac   = (2*r*r_mix - r_mix^2)./r^2;
        alphaf     = alphaf * r_mix_frac;
    else
        r_mix_frac = 1;
        m_s_frac = 1;
    end
    u_eps = alphaf * u * sqrt(rho_B * P.rho_l); % Wrap up the density in u_eps for now
    % ------------------------------
%     toc
    
    % -------- Fragmentation and Grain size output model goes here -------
    % Glass transition limit on fragmentation energy
    hs_sm = (@(T,T0,k) 1./(1 + exp(-2*(3./pI.T_g_rng)*(T-(pI.T_g+pI.T_g_rng/2)))));
    fragPari = pI.fragPari*hs_sm(T_s);
    
%     if T_s > (pI.T_g + pI.T_g_rng)
%         fragPari = pI.fragPari;
%     elseif and(T_s >= pI.T_g, T_s <= (pI.T_g + pI.T_g_rng))
%         fragPari = pI.fragPari * (T_s - pI.T_g)/(pI.T_g_rng);
%     else
%         fragPari = 0*pI.fragPar;
%     end
    fragPar_eff = sum((1-ni).*nsi.*fragPari);                        % Effective fragmentation partitioning, accounting for solids PSD
    
    C_w = heatCapacityP(rho_w,T_w);
    T_i = (T_w * sqrt(pI.rho_m*pI.C_s*pI.k_m) + T_s * sqrt(rho_w * C_w * pI.k_w)) /...
        (sqrt(pI.rho_m*pI.C_s*pI.k_m) + sqrt(rho_w * C_w * pI.k_w));         % Melt interface Temp (Mastin 2007 eq'n 3)
    
    dQ_dSA  = pI.k_m * (T_s - T_i)/pI.r_crit * m_s * m_s_frac / u;         % Heat transfer from particles per surface area [J kg/m^3] 
    
    dm_w_ent = 2 * pi * r * u_eps;  % rho = sqrt(rho*rho_a)? % Water
    
%     [dhw,~] = getDeltaHw(hw,hw0,W0.h_ws,W0.h_vs,pI.C_s,m_s,m_w,dm_w_dz,rho_w,T,W0.Tsat,pI.fragPar);
%     [dhw,~] = getDeltaHw(Hw,hw0,W0.h_ws,W0.h_vs,pI.C_s,m_s*m_s_frac,m_w*m_s_frac,dm_w_dz,rho_w,T_s,W0.Tsat,fragPar_eff);
%     toc
%     tic
    % Estimate max heat transfer in mixing region to bring water/particles
    % to equilibrium
    m_w_mix = dm_w_ent + m_w*m_s_frac;
    Hw_mix = (dm_w_ent*hw0 + m_w*m_s_frac*Hw)/(m_w_mix);
    TeFun = @(x) findTequil(x,pa,m_w_mix,Hw_mix,m_s*m_s_frac,pI.C_s,T_s,fragPar_eff);
    T_w_f = fminbnd(TeFun,pI.Tw0,T_s,opts);
    Hw_f  = enthalpy(density(pa,T_w_f),T_w_f);
    
%     toc
    % Heat transfer from particles limited to max needed to reach equilibrium
%     tic
    dQ_s_dz_max = dQ_dSA * sum(nsi.*pI.SSAi./(1-fragPari));             % Max heat transfer from particles, this time step
    dQ_s_dz = min([dQ_s_dz_max (Hw_f - Hw_mix)*m_w_mix/(1-fragPar_eff)]); %* sum(pI.nsi.*pI.SSAi./(1-fragPari));
    
    m_frag  = dQ_s_dz * fragPar_eff/(1-fragPar_eff) / pI.Ess_out * sum(nsi.*(1-ni));    % fragmentented mass [kg_melt/m]
    nsi_in  = nsi.*fragPari.*pI.SSAi.*(1-ni)./sum(nsi.*fragPari.*pI.SSAi.*(1-ni)); % Distribution of mass removed from coarse fraction
    
    dnsi_dz = m_frag/(m_s*sum((1-ni).*nsi))*(-nsi_in + pI.phi_pdf);      % PSD change
    dm_w_frag = m_frag*(sum((ni).*(nsi)./(1-ni)) - sum((ni).*(nsi+dnsi_dz)./(1-ni)));              % Gas release from fragmentation
    
    % Equilibrate gas release (from pressure change and fragmentation)
%     flargh calc generation of new SA and mass loss per grain size bin
%     dE_s_dz = 1/(1-fragPar) * dQ_s_dz
%     dE_w_dz = dm_w_dz * (gprime_w? * z + hw0) + dQ_s_dz

    gprime_w = (rho_w - rho_wa)./max([rho_w  P.rho_l]) * P.g; % Need persistent
    
%     gprime_s = (rho_s - rho_wa)./max([rho_s  P.rho_l]) * P.g; % Need persistent
    

    % Deltas
    % dm_w_dz, dm_s_dz, dQ_s_dz, dQ_w_dz, dnsi_dz, dpsi_dz
    dm_w_dz = dm_w_ent + dm_w_frag + dm_w_bub;
    dm_s_dz = -dm_w_bub - dm_w_frag;
    dpsi_dz = P.g*(rho_wa - rho_B)*r^2;                             % Momentum
    % dQ_s_dz
%     dQ_w_dz = dm_w_dz * (gprime_w * z + hw0) + dQ_s_dz.*(1-fragPar_eff);
    dE_s_dz = dm_s_dz * (gprime_w * z + H_bw) + dQ_s_dz;
    dE_w_dz = dm_w_ent * (hw0) + (dm_w_bub + dm_w_frag)*(gprime_w * z + H_bw) - dQ_s_dz.*(1-fragPar_eff);
%     toc

    if ~isreal(dm_s_dz) || ~isreal(H_bw) || ~isreal(dQ_s_dz) || ~isreal(gprime_w)
        disp('bollocks')
    end
%     plot(ax(1),z,m_w,'ob')
%     plot(ax(2),dQ_s_dz,'ob')
%     plot(ax(3),T_s,'or')
%     plot(ax(3),T_w,'ob')
    pause(0.02)
%     dE_ss = fragPar/(1-fragPar) * dQ_s_dz / m_s;                    % Specific fragmentation energy, J/kg
%     dA_dz = dE_ss / pI.Es;                                          % Specific surface area, m^2/kg/m

    %   -> A no longer necessary (built into nsi)
    
    % Some sanity checks to test on
%     try
%         assert(~any((nsi+dnsi_dz)<0),'nsi<0')
%     catch
%         disp('Whoa whoa whoa')
%     end
%     toc
%%  Old deltas to be removed
%     [dhw,~] = getDeltaHw(hw,hw0,W0.h_ws,W0.h_vs,pI.C_s,m_s,m_w,dm_w_dz,rho_w,T,W0.Tsat,pI.fragPar);
    % --------
%     flarch check this function ^^ % Compare dhw with (h(Tf)-hw), especially around Tsat
%     dE_ss = pI.fragPar/(1-pI.fragPar)*(dhw) * dm_w_dz/(m_s); % Specific fragmentation energy, J/kg
    
%     gprime = (rho_B - P.rho_l)./max([rho_B P.rho_l]) * P.g; % Reduced grav for bouyancy
%     dE_dz = dm_w_dz * (gprime*z + hw0 - dE_ss);              % Total enthalpy
%     dA_dz = dE_ss / pI.Es;                                % Specific surface area, m^2/kg/m

%%
    % Incremental fragmentation and grain size change COULD happen here
    %   --> Allows incremental gas release from pyroclasts (generally
    %       likely to be negligible at this stage, but could be significant
    %       if introduced all at once later)
    %   --> Allows a parameterization for solids mass loss eventually
    %   --> Requires initial calculation of entrained water, then second
    %   calculation of d_total_water = w_entrain + w_frag_release
    
%     m_w = y(1);
%     m_s = y(2);
%     psi = y(3);
%     Ew  = y(4);
%     Es  = y(5);
%     nsi = y(6:6+numGS-1);
    
    dy_dz    = zeros(5+numGS,1);
    dy_dz(1) = dm_w_dz;
    dy_dz(2) = dm_s_dz;
    dy_dz(3) = dpsi_dz;
    dy_dz(4) = dE_w_dz;
    dy_dz(5) = dE_s_dz;
    dy_dz(6:6+numGS-1) = dnsi_dz;
    
    
    rho_w_prev = rho_w;
%     hw_prev  = hw;
%     Tprev    = T;
%     dm_w_pr  = dm_w_dz;
%     m_w_pr   = m_w;
%     z_pr     = z;
    if nargout>1
        WO.T_s     = T_s;
        WO.T_w     = T_w;
        WO.Tsat    = W0.Tsat;
        WO.rho_w   = rho_w;
        WO.xv      = xv;
        WO.hw      = hw;
    end
%     toc  
        %% --- Old Code for testing fragmentation energy ----
        
%     dE_frag = (h_w_T - h_w0); % dm_w_dz*C_w*(T-Tw0)? % CR*dm_w_dz*Cpm*(T-Tr)*(1-n)?
    % To try a more precise heat transfer estimate for dE_frag
%     Tfun2 = @(x) findT_2(x,pa,T,hw0,pI.C_s,1-pI.fragPar,dm_w_dz,m_s);
%     T2 = fminbnd(Tfun2,Tmin,T);

        
%     if z>65 && z<99
% %         [c_bg,~] = heatCapacityP(P.rho_bg0,T);
% %         c_sb = c_bg*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi); 
%         [dhw_est,Tf_est,x_est] = getDeltaHw(hw,hw0,W0.h_ws,W0.h_vs,pI.C_s,m_s,m_w,dm_w_dz,rho_w,T,W0.Tsat,pI.fragPar);
%         % (T-10)<=W0.Tsat
%         if x_est==0 || x_est==1
%             try 
%                 dhw_act = enthalpy(density(pa,Tf_est),Tf_est) - hw0;
%             catch
%                 disp('eh?')
%             end
%         else
%             dhw_act = enthalpy(densityTX(Tf_est,x_est),Tf_est) - hw0;
%         end
%         dhw_prev = ((1-pI.fragPar)*m_s*pI.C_s*(T-Tprev) + m_w_pr*(hw - hw_prev)) / (-dm_w_pr*(z-z_pr));
% 
%         plot(ax(1),T,hw,'.b')
%         hold(ax(1),'on')
% %         plot(ax(1),T,hw0+dhw_act,'ks')
%         plot(ax(1),T,hw0+dhw_est,'or')
%         plot(ax(1),T,hw0+dhw_prev,'om')
%         plot(ax(2),T,xv,'.b')
% %         plot(ax(2),dm_w_pr,dhw_prev,'or')
%         hold(ax(2),'on')
% %         plot(ax(2),dm_w_dz,dhw_est,'.b')
%         plot(ax(2),T,x_est,'or')
% %         plot(ax(2),T,T-Tprev,'or')
% %         plot(ax(3),T,dhw_act,'.b')
%         plot(ax(3),T,dhw_prev/dhw_est,'.b')
%         hold(ax(3),'on')
% %         plot(ax(3),T,dhw_est,'or')
% %         plot(ax(3),T,dhw_prev,'ks')
%         pause(0.1)
%     elseif z>99
%         ylabel(ax(1),'h_w')
%         legend(ax(1),{'Current h_{w,i}','hw_0+dh_w','i-1'})
%         ylabel(ax(2),'x')
%         legend(ax(2),{'xv','xv_{est}'})
%         ylabel(ax(3),'dh_{w,i-1}/dh_{w,est}')   
%     end

end

function [value,isterminal,direction] = stopEntrain(y,m_s)
    % u = 0, T ~ Tw0?
    m_w = y(2);
    psi = y(3);
    m   = m_w + m_s;
    u   = psi/m;

    value      = u;
    isterminal = 1;
    direction  = 0;

end

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
        Tf = ((1-f)*m_s*c_s*T + m_w*cpv*T - dm_w*(hw - cpv*T - hw0)) / (cpv*(m_w+dm_w) + (1-f)*m_s*c_s);
        dhw = hw + cpw*(Tf - T) - hw0;
        x = 0;
    end
    
end

function SSA = SSAfromD(ri,rhoi,D,A)
% Calculate specific surface area, m^2/kg, from power-law grain size distribution
%   ri   = vector of effective particle radii
%   rhoi = density of particles at each radii
%   D    = power law exponent
%   A    = roughness scale parameter

    SSA = (3*A*sum(ri.^(2-D))./sum(rhoi.*ri.^(3-D)));

end