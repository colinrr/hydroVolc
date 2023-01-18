function [pI2,wO] = MWIv1(cI,cO,pI,dO)
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
%     global xv
    P.xv_crit     = 0.005; % Dryness fraction below which to stop run and adjust tolerances
    
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
    
    % Enthalpy
%     C_sb0   = W0.C_p.*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi); % Bulk pyroclast heat capacity
%     H0      = n0*W0.h + (1-n0)*(C_sb0)*(pI.T0 - P.Tr);
    % Alt Enthalpy - more accurate but small diff
    h_s     = enthalpy(P.rho_bg0,pI.T0)*sum(pI.ni.*pI.nsi) ...
                + pI.C_s*sum((1-pI.ni).*pI.nsi)*(pI.T0 - P.Tr);
    H0      = n0*W0.h + (1-n0)*h_s;
    
    gprime  = (pI.rho_B0 - P.rho_l)./max([pI.rho_B0 P.rho_l]) * P.g;
    E0      = m_0 * (H0 + gprime*z0 + 0.5*pI.u_0^2); % Should this be g'* z' to give positive bouyancy in water?
    
    SSAfromD = @(ri,rhoi,D,A) (3*A*sum(ri.^(2-D))./sum(rhoi.*ri.^(3-D)));
    SSA0    = SSAfromD(pI.Rgsd,pI.rhoi,pI.D,pI.saScale);
    % A0 - Total particle surface area. SSA should be enough?

    % u0 = pI.u_0;
    
    P.Ri0       = P.g*((P.rho_l-pI.rho_B0)/P.rho_l)*pI.r_0/(pI.u_0^2);
    P.De        = 2*pI.r_0;
    P.Lme       = pi^0.25*pI.r_0*abs(P.Ri0)^(-0.5);
    
    % Jet entrance scaling
%     gprime = P.g*(rho_wa-rho_B)/max([rho_wa rho_B]);
    Mo =  pi*(pI.r_0).^2.*pI.u_0.^2;
%     Fo =  pi*gprime*(pI.r_0).^2*pI.u_0;
%     Lj =  P.LjC * Mo^(3/4) * Fo^(-1/2);
    P.Lj        = Mo^(1/4)*(pI.u_0/pI.r_0).^(1/2);
    
%% SOLVE ODE

    % ============ intitial run ===========
    options = odeset('RelTol',1e-6,'AbsTol',1e-3,'Events',@(z,y) stopEntrain1(y,m_s0,P.xv_crit));
    zSpan = [pI.Ld cI.Zw]; % Z reference frame: Ld = 0? vent = 0?
    IC    = [m_w0 psi0 E0 SSA0];
%     tic
    odeSolution = ode15s(@(z,y) waterEntrain(z,y,cI,pI,P),zSpan,IC,options);
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
    m_s = m_s0; % COULD include a sedimentation parameterization in future?
    wO.psi = Y(:,2);
    wO.E   = Y(:,3);
    wO.A   = Y(:,4);

    % ---- Post-process -----
    pI2 = pI;
    
    wO.z = Z;
    wO.m = m_s + m_w;
    wO.u = wO.psi./wO.m;
    wO.H = wO.E./wO.m;
    wO.P = P.rho_l.*P.g.*(cI.Zw - Z) + P.Patmo;
    
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
        wO.xv(zi)    = W.xv;
        wO.hw(zi)    = W.hw;
        wO.alpha(zi) = W.alpha;
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
    D_all   = (2.8:.02:4);
    SSA_all = SSAfromD(pI.Rgsd,pI.rhoi,D_all,pI.saScale);
    wO.D    = interp1(SSA_all,D_all,wO.A(end),'spline');
    pI2.D   = wO.D;
    [pI2.nsi,pI2.ni,pI2.Rgsd,pI2.rhoi,pI2.porR] = getPSD(cI.rho_melt,pI.phi0, cO.pg(end), wO.T(end), wO.xv(end), wO.D,pI.phiSz_min,pI.phiSz_max);

    
    % Update plume params    
    pI2.n_0     = (m_w(end) + m_s.*(sum(pI.ni.*pI.nsi) - sum(pI2.ni.*pI2.nsi))) ./ (m_w(end) + m_s);
    pI2.T0      = wO.T(end);
    pI2.vh0     = cI.vh0 + wO.z(end);
    pI2.u_0     = wO.u(end);
    pI2.rho_B0  = wO.rho_B(1);
    pI2.r_0     = ((m_s + m_w(end))./(pi*wO.rho_B(end)*pI.u_0)).^(1/2); % PROBLEM WITH PLUME RADIAL EXPANSION THAT IS FAR TOO RAPID HERE
    pI2.xv_0    = wO.xv(end);
    
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
    h_s_sat = h_bg*n_gs + C_s*(1-n_gs)*(T - Tr); 
    
    H_est = (1-n).*h_s_sat + n.*hw;
    dT = abs(H - H_est);
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

persistent rho_prev % zchk % hw_prev Tprev dm_w_pr m_w_pr ax z_pr
global xv
if isempty(rho_prev)
    rho_prev = pI.rho_B0;
%     zchk = false;
%     tic
%     figure(1)
%     ax(1)=subplot(3,1,1);
%     ax(2)=subplot(3,1,2);
%     ax(3)=subplot(3,1,3);

end
%     zvals = [70:10:170];
%     zchk = false;
%     if any(abs(zvals-z)<1) && ~zchk
%         tic
%     end


%     z = y(1);
    m_w = y(1);
    m_s = cI.Q*(1-pI.n_0);
    psi = y(2);
    E   = y(3);
    A   = y(4);
    
    m = m_s + m_w;
    n = m_w/m;
    u = psi/m;
    gprime = (rho_prev - P.rho_l)./max([rho_prev P.rho_l]) * P.g;
    H = E/m - u^2/2 - gprime*z;
    pa = P.rho_l.*P.g.*(cI.Zw - z) + P.Patmo;  % Ambient pressure

    % ----- Get water properties and mixture Temperature ------
    W0 = waterProps(pa,pI.T0); 
%     C_sb = W0.C_p.*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi); % Bulk pyroclast heat capacity
    Tmax = H/((1-n)*pI.C_s) + P.Tr; % Use C_b here?
%     Tmax = min([H/((1-n)*pI.C_s) + P.Tr pI.T0]); % Use C_b here?

%     C_sb   = W0.C_p.*sum(pI.ni.*pI.nsi) + pI.C_s*sum((1-pI.ni).*pI.nsi)
%     Hvap = (1-n)*C_sb*(W0.Tsat - P.Tr) + n*W0.h_vs;
%     Hliq = (1-n)*C_sb*(W0.Tsat - P.Tr) + n*W0.h_ws;

    % Pyroclast enthalpy
    % --> Bubble gas is minor component, assume isochoric process for now. 
    %     Room here for interesting work on pyroclast water ingestion and 
    %     related effects
%     pb = pressure(P.rho_bg0,W0.Tsat);
    [h_bgS,~] = enthalpy(P.rho_bg0,W0.Tsat);
    h_s_sat = h_bgS*sum(pI.ni.*pI.nsi) ...
            + pI.C_s*sum((1-pI.ni).*pI.nsi)*(W0.Tsat - P.Tr); 
        
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
      % !!!NOTE!!!: Possible that the linear transition to xv=0 near
      % h_ws is causing the integrator to fail. Possibly try a gaussian
      % smooth here?
    % ---------------------------------------------------------

    
    hw      = enthalpy(rho_w,T);                    % Plume water enthalpy
    rho_wa  = density(pa,pI.Tw0);
    hw0     = enthalpy(rho_wa,pI.Tw0);              % Ambient water enthalpy
    
    rho_s   = (sum(pI.nsi./pI.rhoi))^-1;            % Solids bulk density (assuming no change in D)
    rho_B   = (n/rho_w + (1-n)/rho_s)^-1;           % Plume bulk density
    r       = (m/(rho_B*u))^(1/2);                  % Plume radius
    
    % ---- Entrainment model ------
    Ri      = P.g*((P.rho_l-rho_B)/P.rho_l)*r/(u^2); % Local Richardson #
    alphaf  = alpha_cara(Ri,z,P.De,P.Lme,r);
    if alphaf<0.05
        alphaf=0.05;
    elseif alphaf>0.17
        alphaf=0.17;
    end
    
    % Simple scale for mixing layer development over a jet entrance length
    if and( pI.useJetEntranceLength , (z-pI.Ld)<P.Lj)
        alphaf = alphaf * (z-pI.Ld)./P.Lj; 
    end
    u_eps = alphaf * u * sqrt(rho_B * P.rho_l); % Wrap up the density in u_eps for now
    % ------------------------------
    
    % Deltas

    dm_w_dz = 2 * pi * r * u_eps;  % rho = sqrt(rho*rho_a)? % Water
    dpsi_dz = P.g*(P.rho_l - rho_B)*r^2;                    % Momentum
    
    
%     dE_frag = (h_w_T - h_w0); % dm_w_dz*C_w*(T-Tw0)? % CR*dm_w_dz*Cpm*(T-Tr)*(1-n)?
    % To try a more precise heat transfer estimate for dE_frag
%     Tfun2 = @(x) findT_2(x,pa,T,hw0,pI.C_s,1-pI.fragPar,dm_w_dz,m_s);
%     T2 = fminbnd(Tfun2,Tmin,T);

    %% --- Code for testing fragmentation energy ----
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
%         plot(ax(2),T,X,'.b')
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
%         legend(ax(2),{'X','X_{est}'})
%         ylabel(ax(3),'dh_{w,i-1}/dh_{w,est}')   
%     end
%%
    [dhw,~] = getDeltaHw(hw,hw0,W0.h_ws,W0.h_vs,pI.C_s,m_s,m_w,dm_w_dz,rho_w,T,W0.Tsat,pI.fragPar);
    % --------
    
    % Glass transition limit on fragmentation energy
    if T > (pI.T_g + pI.T_g_rng)
        fragPar = pI.fragPar;
    elseif and(T >= pI.T_g, T <= (pI.T_g + pI.T_g_rng))
        fragPar = pI.fragPar * (T - pI.T_g)/(pI.T_g_rng);
    else
        fragPar = 0;
    end
    dE_ss = fragPar/(1-fragPar)*(dhw) * dm_w_dz/(m_s); % Specific fragmentation energy, J/kg_solids
    
    gprime = (rho_B - P.rho_l)./max([rho_B P.rho_l]) * P.g; % Reduced grav for bouyancy
    dE_dz = dm_w_dz * (gprime*z + hw0) - m_s*dE_ss;              % Total enthalpy change
    dA_dz = dE_ss / pI.Es;                                % Specific surface area, m^2/kg_solids/m

    % Incremental fragmentation and grain size change COULD happen here
    %   --> Allows incremental gas release from pyroclasts (generally
    %       likely to be negligible at this stage, but could be significant
    %       if introduced all at once later)
    %   --> Allows a parameterization for solids mass loss later
    %   --> Requires initial calculation of entrained water, then second
    %   calculation of d_total_water = w_entrain + w_frag_release
    
    dy_dz    = zeros(4,1);
    dy_dz(1) = dm_w_dz;
    dy_dz(2) = dpsi_dz;
    dy_dz(3) = dE_dz;
    dy_dz(4) = dA_dz;
    
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

function [value,isterminal,direction] = stopEntrain1(y,m_s,xv_crit)
    % Stop integration of velocity reaches zero, or water dryness fraction
    % gets close to 0

    global xv
    
    % u = 0, T ~ Tw0?
    m_w = y(2);
    psi = y(3);
    m   = m_w + m_s;
    u   = psi/m;

    value      = [u xv-xv_crit];
    isterminal = [1 1];
    direction  = [0 -1];

end

function [value,isterminal,direction] = stopEntrain2(y,m_s)
    % Stop integration of velocity reaches zero only

    
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
        Tf = ((1-f)*m_s*c_s*T + m_w*cpw*T - dm_w*(hw - cpw*T - hw0)) / (cpw*(m_w+dm_w) + (1-f)*m_s*c_s);
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