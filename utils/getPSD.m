function [nsi,ni,r,rho_s,porosity] = getPSD(rho_m,porosity0, Pg, T, X, D,phi_min,phi_max)
% Generate grain size distribution given power law exponent and size range,
% base on Girault 2014 and Kaminski and Jaupart 1998.
% IN:
%   rho_m       = melt density (kg/m^3)
%   porosity0   = pyroclast porosity (large particles)
%   Pg          = gas pressure in vesicles (Pa)
%   T           = gas temperature in vesicles (K)
%   X           = dryness fraction (default = 1) for water at saturation
% OPTIONAL IN:
%   D           = power law exponent. Default 3.0
%   phi_min     = minimum phi size (max grain size)
%   phi_max     = max phi size (min grain size)
% OUT:
%   nsi          = particle mass fraction at each grain size (including
%                 bubble gas). Sum(nsi) = 1.
%   ni          = mass fraction of bubbles at each grain size, relative to
%                 particle mass fraction.
%   r           = particle radii (m)
%   rho_s       = bulk pyroclast density as function of grain size
%   porosity    = output porosity as a function of grain size
%
% C Rowell Mar 2021

narginchk(4,8)
if nargin<8
    phi_max = 10;
end
if nargin<7
    phi_min = -9;
end
if nargin<6
    D = 3.0;
end
if nargin<5
    X = 1;
end

if isempty(phi_max)
    phi_max = 10;
end
if isempty(phi_min)
    phi_min = -9;
end
if isempty(D)
    D = 3.1;
end

    Rv   = 461; % Gas constant for water vapour
%%
    N0  = 1;
    % rho_m = 2400;

    rho_b = rho_m;

    % Size cutoffs for porosity - Kaminski and Jaupart 1998
    rbig = .01;  % radius above which effective porosity is flat. 1 cm, KJ '98
    rc = 1e-4;  % critical radius for no effective porosity. 100 um, KJ '98


    phi = (phi_min:1:phi_max)';    % Phi vector
    r   = (2.^-phi)/2000;         % Radius (m)
    nr  = length(r);

    % Set up porosity and density as f(r)
    rsmall = r<rc;
    rmid   = and(r>=rc,r<rbig);
    porosity = porosity0*ones(nr,1);
    porosity(rmid) = porosity(rmid).*(1-rc./r(rmid));
    porosity(rsmall) = zeros(sum(rsmall),1);

    % Particle density as a function of radius
%     rho_s = (1-porosity)*rho_m + (porosity)*Pg/(Rv*T);
%     rho_s = (1-porosity)*rho_m + (porosity)*EoS_H2O_2(Pg,T);
    if X<1
        rho_w = densityTX(T,X);
    else
        rho_w = density(Pg,T);
    end
    rho_s = (1-porosity)*rho_m + (porosity)*rho_w;


    % Get power-law GSD
%     Nr = flipud(2.^(log2(N0) - D.*phi));    % Cumulative number (R>=r)
%     Ni = [Nr(1); Nr(2:end) - Nr(1:end-1)];  % Number per bin
    Ni = 2.^(log2(N0) + D.*phi);

    nsi = Ni.*r.^3.*rho_s./sum(Ni.*r.^3.*rho_s); % Mass fractions 
    
    % Likely also want a gas mass fraction/ratio at each grain
    % size...assuming none in melt for the moment
%     C = porosity./(1-porosity).*Pg./(Rv*T*rho_m);
    C = porosity./(1-porosity).*rho_w/rho_m;
    ni = C./(1+C);

    % semilogx(r,mi)
%     plot(phi,mi)
end