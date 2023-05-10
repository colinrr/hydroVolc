function [N2,yz] = getBruntVaisala(Z,P,T,RH,z,zRef,Nmode)
% Given atmo profile P(Z), T(Z), RH(Z), calculate squared Brunt Vaisala 
% frequency at heights z
% Z -> profile heights (meters)
% P -> profile pressures (Pascals)
% T -> profile Temperature (Kelvin)
% RH-> profile relative humidity (fraction, not percent)
%
% z -> Optional altitude values for output N [default is all Z]
%       []     : get N.^2(Z) for whole profile Z
%       z1     : scalar z to get a single local gradient value at height z
%       [z1 z2...zn]: get average N.^2 between each height (so long as n<length(Z))
%
% zRef -> Height for reference density or pressure
%       []: Default uses 1000 hPa reference P (for potential T mode) 
%           or local density rho(z) for density gradient mode
%       scalar: calculates N.^2(z), using P(z=zRef) or rho(z=zRef)
%       'tp'   :    uses density at height of tropopause as the reference
%
% Nmode -> Calculation method
%           'theta': Use potential Temperature (default)
%           'drho':  Linear density gradient approximationa
%           'gas':   Use gas constant formulation (no reference P or rho)
%
% OUTPUT:
%   N2:     Squared BV frequency (1/s.^2)
%   yz:     profile as function of z. 
%             ->  In 'theta' mode, this is potential temperature theta (K)
%             ->  In 'drho' and 'gas' modes, this is bulk atmospheric
%                 density, rho_aB (kg/m^3)

narginchk(4,7)
assert(isvector(Z))
assert(isvector(P))
assert(isvector(T))
assert(isvector(RH))
assert(length(P)==length(Z))
assert(length(T)==length(Z))
assert(length(RH)==length(Z))

if nargin<5
    z = [];
end
if nargin<6
    zRef = []; %max([min(Z) 0]);
end
if nargin<7
    Nmode = 'theta';
end
if ~any([strcmp(Nmode,'theta') strcmp(Nmode,'drho') strcmp(Nmode,'gas')])
    Nmode = 'theta';
    warning('Brunt-Vaisalla: default to potential temperature formula.')
end

if isempty(z)
    zMode = 'all';
    z = Z;
elseif isvector(z) && isnumeric(z) && length(z)<length(Z)
    zMode = 'avg';
else
    error('zMode not recognized, something broke')
end

if isempty(zRef)
    zRefMode = 'local';
elseif isscalar(zRef) && isnumeric(zRef)
    zRefMode = 'fixed';
    Pref = [];
elseif ischar(zRef) && strcmp(zRef,'tp')
    zRef = findTPheight(Z/1e3,T);
    zRefMode = 'fixed';
    Pref = [];
else
    error('zRefMode not recognized.')
end

    g       = 9.81;  % gravitational acceleration (m/s^2)
    R_v     = 461;   % gas constant of volcanic gas (water) (J/kg/K) 
    R_d     = 287;   % gas constant of dry air (J/kg/K) 
    Cpd     = 1005;  % Air constant pressure heat capacity (J/kg/K)
    Cvd     = 718;   % Air constant volume heat capacity   (J/kg/K)
    eps       = R_d/R_v;

    % Get atmospheric density profile
    es = 100.*6.112*exp(17.67*(T-273.15)./(T +243.5-273.15)); % saturation vapor pressure for water
    dp = max([P-es 2*es],[],2); % This is janky at such low pressures
    w_s = 1./eps * es./dp; %(P-es); % w_a = r_h*w_s, mass mixing ratio of water vapor to dry air at saturation
    w_a = RH.*w_s; % mass mixing ratio of water vapor to dry air

    Cpv = 1860; % Specific heat, water vapor at 273 K (assuming constant for now)
    Cvv = Cpv - R_v;
    Cpa = (Cpd + w_a*Cpv)./(1+w_a);  % Atmo bulk specific heat capacity constant P
    Cva = (Cvd + w_a*Cvv)./(1+w_a);  % Atmo bulk specific heat capacity constant V

    if strcmp(zRefMode,'tp')
        zRef = findTPheight(Z/1e3,T);
    end
%     gamma = Cpa./Cva;
%     theta = T.*(P(1)./P).^((Cpa-Cva)./Cpa); % Potential Temp.

    

    switch Nmode
        case 'theta'
            % Get reference P    
            switch zRefMode
                case 'local'
                    Pref = 1000; % hPa
                case 'fixed'
                    P = 10^interp1(Z,log10(P),real(zRef),'pchip','extrap');
%                 case 'tp'              
            end
            theta = T.*(Pref./P).^((Cpa-Cva)./Cpa); % Potential Temp.
            N2 = g./theta .* gradient(theta,Z);
            
            yz = theta;

        case 'drho'
            rho_aB    = P./(R_v.*T).*(1+w_a)./(w_a+eps);  % Atmospheric gas (air+vapour) bulk density
            drho_dZ   = gradient(rho_aB,Z);           % Density gradient
            switch zRefMode
                case 'local'
                    rhoRef = rho_aB;
                case 'fixed'
                    rhoRef = 10.^interp1(Z,rho_aB,zRef,'spline','extrap');
            end
            N2 = -g./rhoRef .* drho_dZ;
            
            yz = rho_aB;
            
        case 'gas'
            rho_aB    = P./(R_v.*T).*(1+w_a)./(w_a+eps);  % Atmospheric gas (air+vapour) bulk density
            drho_dZ   = gradient(rho_aB,Z);           % Density gradient
            gamma = Cpa./Cva;
            N2 = g.*(1./gamma .* gradient(log(P),Z) - gradient(log(rho_aB),Z)); % Alternative formulation using gas constants, local only

            yz = rho_aB;

    end
    
    switch zMode
        case 'avg'
            N2avg = zeros(length(z)-1,1);
            for zi=1:(length(z)-1)
                zii = and(Z>=z(zi),Z<=z(zi+1));
                if sum(zii)==0
                    error('Brunt Vaisalla averaging interval is not within bounds.');
                elseif sum(zii)<100
                    zt = linspace(z(zi),z(zi+1),100);
                    Nt = interp1(Z,N2,zt,'pchip');
                else
                    zt = Z(zii);
                    Nt = N2(zii);
                end
                N2avg(zi) = trapz(Z(zii),N2(zii))./diff(z(zi:zi+1)); % integral average
            end
            N2 = N2avg;
    end
    
