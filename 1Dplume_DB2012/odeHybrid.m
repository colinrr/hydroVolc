function dy_ds = odeHybridtest(y, ...
    g, ...
    C_d, ...
    C_v, ...
    C_l, ...
    C_s, ...
    R_d, ...
    R_v, ...
    L,   ...
    model,      ...
    nexp,      ...
    alpha,      ...
    beta,       ...
    baratio,       ...
    omega,      ...
    mu_a,   ...
    ri,     ...
    rhoi,   ...  %  m_s0,       ...
    ni,   ...
    De,     ...
    Lme,     ...
    u_0,       ...
    vent_height,       ...
    Meteo_Humidity,        ...
    Meteo_Height, ...
    Meteo_Temperature, ...
    Meteo_Wind, ...
    Meteo_Pressure)




if abs(real(y(2))./imag(y(2)))<10 
'complex value calculated'
a=thisvariabledoesnotexist;
end
    numGS = length(ri);
    
    z         = y(2);
    m_d       = y(3);
    m_v       = y(4);
    m_l       = y(5);
%     m_s       = m_s0;
    m_si      = y(6:6+numGS-1);
    psi       = y(6+numGS);
    angle     = y(7+numGS);
    Q         = y(8+numGS);
    % P         = y(9); uncomment if hydrostatic approx used
    if real(z)<=max(Meteo_Height)
        Va = interp1(Meteo_Height,Meteo_Wind,real(z),'pchip','extrap');
    else
        Va = interp1(Meteo_Height,Meteo_Wind,real(z),'pchip',Meteo_Wind(end));
    end
    
    if Va<0
        Va=0;
    end
    if Va>200
        Va=200;
    end
    
    wse=Va/u_0;
    
    
    
    if real(z)<=max(Meteo_Height)
        theta_a = interp1(Meteo_Height,Meteo_Temperature,real(z),'pchip','extrap');
    else
        theta_a = interp1([Meteo_Height;50000-vent_height;85000-vent_height],[Meteo_Temperature;273;173],real(z),'pchip','extrap');
    end
    
    
    
    % Ambient humidity and pressure -----
    hur_a=interp1(Meteo_Height,Meteo_Humidity,real(z),'pchip','extrap');
    if real(z)>max(Meteo_Height)
        hur_a=0;
    else
        hur_a=interp1(Meteo_Height,Meteo_Humidity,real(z),'pchip','extrap');
    end
    
    if hur_a<0
        hur_a=0;
    elseif hur_a>1
        hur_a=1;
    end
    
    P = 10^interp1(Meteo_Height,log10(Meteo_Pressure),real(z),'pchip','extrap');
    % -----
    
    eps       = R_d/R_v;
    
    
    % Atmospheric vapor
    es = 100.*6.112*exp(17.67*(theta_a-273.15)./(theta_a +243.5-273.15)); % saturation vapor pressure for water
    dp = max([P-es 2*es]);
    w_s = 1/eps * es/dp; %(P-es); % w_a = r_h*w_s, mass mixing ratio of water vapor to dry air at satuartion
%     w_s = 1/eps * es/(P-es); % w_a = r_h*w_s, mass mixing ratio of water vapor to dry air at satuartion
    w_a = hur_a*w_s; % mass mixing ratio of water vapor to dry air
    
    
    % set up ODE's
    
    rho_aB    = P/(R_v*theta_a)*(1+w_a)/(w_a+eps);  % Atmospheric gas (air+vapour) bulk density
    rhophi_av = P/(R_v*theta_a)*w_a/(w_a+eps);      % Atmo vapour volume fraction
    rhophi_ad = P/(R_v*theta_a)*1/(w_a+eps);        % Atmo dry air volume fraction
    C_aB      = (C_d + w_a*C_v)/(1+w_a);            % Atmo bulk specific heat capacity
    
    % -- For GSD, CR Mar 2021 ---
    m_s       = sum(m_si);                          % Total pyroclast mass
    nsi        = m_si./m_s;                          % Pyroclast mass fractions by grain size
    % --
    
    m         = m_d + m_v + m_l + m_s;              % Total mass flux
    u         = psi/m;                              % Velocity from momentum/mass
    % display(m)
    C_sb      = C_v.*sum(ni.*nsi) + C_s*sum((1-ni).*nsi); % Bulk pyroclast heat capacity
    C_B       = (m_d*C_d + m_v*C_v + m_l*C_l + m_s*C_sb)/m; % Bulk mixture heat capacity
    % display(C_B)
    theta     = (1/C_B)*(Q/m);                      % Plume bulk temperature
    rho_B     = (P/R_v/theta)*m/(m_v+eps*m_d);      % Plume bulk density
    r         = (m/(rho_B*u))^(1/2);                % Plume radius
    rhophi_v  = (P/R_v/theta)*m_v/(m_v+eps*m_d);    % Bulk density of vapor fraction
    
    
    % entrainment assumption
    % if rho_B > rho_aB
    %     u_eps    = 1/16*(rho_B/rho_aB)^0.5*abs(u-Va*cos(angle)) + beta*abs(Va*sin(angle));
    % %     u_eps    = 0.05*(rho_B/rho_aB)^0.5*abs(u-Va*cos(angle)) + beta*abs(Va*sin(angle));
    % else
    Ri=g*((rho_aB-rho_B)/rho_aB)*r/(u^2);
    if model==1
        alphaf=alpha;betaf=beta;
    elseif model==2
        alphaf=alpha_cara(Ri*sin(angle),z,De,Lme,r);betaf=beta;
    elseif model==3
        alphaf=alpha_cara(Ri*sin(angle),z,De,Lme,r);betaf=beta_folch(Ri,wse);
    elseif model==4
        alphaf=alpha_cara(Ri*sin(angle),z,De,Lme,r);betaf=baratio*alphaf;
    else
        'unknown model number'
        tototo=tatata
    end
    
    % figure(1)
    % plot(alphaf,z,'rx')
    % hold on
    if alphaf<0.05
        alphaf=0.05;
    end
    if alphaf>0.17
        alphaf=0.17;
    end
    if betaf<0.1
        betaf=0.1;
    end
    if betaf>1
        betaf=1;
    end
    
    % Entrainment velocity
    u_eps    = ((alphaf*abs(u-Va*cos(angle)))^nexp + (betaf*abs(Va*sin(angle)))^nexp)^(1/nexp);
    % end

    % --- GSD, CR Mar 2021 ---
    zeta    = 0.27; %2*((1+6/5*alphaf).^2 - 1)./((1+6/5*alphaf).^2 + 1);        % Fallout probability
    uf_all  = [((6.2.*ri.*g.*(rhoi-rho_B)/rho_aB).^(1/2))',                     % Fallout velocity at all Re_p
               (ri.*(4.*g^2*(rhoi-rho_B).^2/(225*mu_a*rho_aB)).^(1/3))',
               (4.*ri.^2.*g.*(rhoi-rho_B)/(18*mu_a))']';
           
    Re_p    = 2*ri.*uf_all*rho_aB/mu_a;                                         % Particle Reynold's numbers
    ucheck  = [Re_p(:,1)>=500 and(Re_p(:,2)>=0.4,Re_p(:,2)<500) Re_p(:,3)<0.4]; % Check u_fall and Re overlap
    ucheck  = [ucheck(:,1) and(~ucheck(:,1),ucheck(:,2)) and(~or(ucheck(:,1),ucheck(:,2)),ucheck(:,3) )]; % Resolve ambiguities towards high Re_p
    uf      = uf_all(ucheck);

%         assert(length(uf)==size(uf_all,1),'Ambiguous calculation for settling velocity.')
    
    dm_si_ds  = -zeta*m_si.*uf./(r*u);                                          % Rate of particle mass loss
    % -------
    
    dx_ds     = cos(angle);
    dz_ds     = sin(angle);
    dm_d_ds   = 2*u_eps*r*rhophi_ad;
    dm_v_ds   = 2*u_eps*r*rhophi_av - omega*rhophi_v*r^2;
    dm_l_ds   = omega*rhophi_v*r^2;
    dm_f_ds   = dm_d_ds + dm_v_ds + dm_l_ds;
    dm_ds     = dm_d_ds + dm_v_ds + dm_l_ds + sum(dm_si_ds);
    dpsi_ds   = g*(rho_aB - rho_B)*r^2*sin(angle) + Va*cos(angle)*dm_ds + u.*sum(dm_si_ds);
    dangle_ds = (1/psi)*(g*(rho_aB - rho_B)*r^2*cos(angle) - Va*sin(angle)*dm_ds);
    dQ_ds     = C_aB*theta_a*dm_f_ds - rho_B*u*r^2*g*sin(angle) + L*dm_l_ds + C_sb*theta.*sum(dm_si_ds);

    
    dy_ds=zeros(8,1); 
    dy_ds(1) = dx_ds;
    dy_ds(2) = dz_ds;
    dy_ds(3) = dm_d_ds;
    dy_ds(4) = dm_v_ds;
    dy_ds(5) = dm_l_ds;
    dy_ds(6:6+numGS-1)  = dm_si_ds;
    dy_ds(6+numGS)      = dpsi_ds;
    dy_ds(7+numGS)      = dangle_ds;
    dy_ds(8+numGS)      = dQ_ds;

% dy_ds(9) = dP_ds;