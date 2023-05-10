function [Ri,Ri0] = getRichardsonProfile(dat,z, RiMode)
% Get Richardson numbers for plume model run.
% dat =   Model run output struct (scalar)
% z   =   Optional height to for calculation (default is all z)
% RiMode = 'simple': use simple gprime * r / u.^2 approximation
%          'gradient': gradient Richardson number [default]
%          'BV': Use Brunt-Vaisalla formulation with potential temperature
%
% OUT:
%  Ri  = Ri(z) local Richardson number with height
%  Ri0 = Source Richardson number at plume starting height (vent or water surface
%        level)

narginchk(1,3)
if nargin<3
    RiMode = 'gradient';
end
if nargin<2 || isempty(z)
    z = dat.pO.z;
    interpz = false;
else
    interpz = true;
end
if ~any([strcmp(RiMode,'simple') strcmp(RiMode,'gradient') strcmp(RiMode,'BV')])
    RiMode = 'gradient';
    warning('Richardson: default to gradient formula.')
end

g = 9.81;

switch RiMode
    % Simple formulation, gprime * r / u.^2, typical for source experiments, but not great when the flow is not Boussinesq
    case 'simple' 
        gprime = g.* (dat.pO.rho_B - dat.pO.atmo.rho)./max([dat.pO.rho_B dat.pO.atmo.rho],[],2);
        Ri = gprime .* dat.pO.r ./ dat.pO.u.^2;
    
    % More general gradient Richardson number
    case 'gradient' 
        Ri = g./ dat.pO.rho_B .* gradient(dat.pO.rho_B,dat.pO.z) ./ gradient(dat.pO.u,dat.pO.z);

    % Brunt-Vaisalla formulation (useful when using potential temperature 
    % to get atmospheric conditions, but does not account for plume local density)
    case 'BV'
        atmo = interpAtmoArray(dat.pI.atmo,dat.pO.z);
        N2 = getBruntVaisala(dat.pO.z,atmo(:,1),atmo(:,3),...
    atmo(:,5)/100);
        Ri = N2 ./ gradient(dat.pO.u,dat.pO.z);
        
end

Ri0 = Ri(1);
if interpz
    Ri = interp1(dat.pO.z,Ri,z,'pchip','extrap');
end

end