function aP = getAtmoProps(atmo,z)
    % z = m a.s.l. vector
    % atmo = amto profile
    
    if(ischar(atmo))
        load(atmo)
        clear atmo
        if or(istable(atmprofile),isstruct(atmprofile))
            atmo.z  = atmprofile.Altitude;
            atmo.P  = atmprofile.Pressure;
            atmo.T  = atmprofile.Temperature;
            atmo.rh = atmprofile.relativeHumidity/100;
        end
    elseif ismatrix(atmo) && size(atmo,2)==5
        A = atmo; clear atmo
        atmo.z = A(:,2);
        atmo.P = A(:,1);
        atmo.T = A(:,3);
        atmo.rh = A(:,5)/100;
    end
    
    R_d        = 287;   % gas constant of dry air (J/kg/K) 
    R_v        = 461;   % gas constant of volcanic gas (water) (J/kg/K) 
    eps        = R_d/R_v;
    g          = 9.81;
    
    aP.z    = z;
    aP.P    = 10.^interp1(atmo.z,log10(atmo.P),real(z),'pchip','extrap');
    if real(z)<=max(atmo.z)
        aP.T = interp1(atmo.z,atmo.T,real(z),'pchip','extrap');
    else
        aP.T = interp1([atmo.z;50000;85000],[atmo.T;273;173],real(z),'pchip','extrap');
    end
    aP.rh   = interp1(atmo.z,atmo.rh,real(z),'pchip','extrap');
    aP.rh(aP.rh<0) = 0;
    
    aP.es_0   = 100.*6.112.*exp(17.67.*(aP.T-273.15)./(aP.T +243.5-273.15)); % saturation vapor pressure for water
    dp = max([aP.P-aP.es_0  2*aP.es_0],[],2); % Avoid a div0 when abs. pressure becomes very small
    aP.w_s    = 1./eps .* aP.es_0 ./ dp; % w_a = r_h*w_s, mass mixing ratio of water vapor to dry air at saturation
    aP.w_a    = aP.rh.*aP.w_s; % mass mixing ratio of water vapor to dry air
    aP.rho_aB = aP.P./(R_v.*aP.T).*(1+aP.w_a)./(aP.w_a+eps);  % Atmospheric gas (air+vapour) bulk density

    aP.dpdz   = gradient(aP.P,z);
    aP.drhodz = gradient(aP.rho_aB,z);
    
    % Calculate Dew-point Temperature
    qz      =   aP.rh .* 0.622 .* aP.es_0 ./ (aP.P - aP.es_0);
    chi     = log(aP.P./100 .* qz ./ (6.112.*(0.622+qz)));
    aP.Tdew = 243.5 .* chi ./ (17.67-chi);
    
    % Local Brunt-Vaisalla
%     N = (-g./aP.rho_aB .* aP.drhodz).^(1/2); % uses a local reference density?
    
end