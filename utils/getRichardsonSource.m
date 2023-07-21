function Ri0 = getRichardsonSource(qA, RiMode)
% Get source Richardson numbers for plume model sweep.
% qA  =   Model sweep summary (see getHydroVolcSweepSummary)
% z   =   Optional height to for calculation (default is all z)
% RiMode = 'simple': use simple gprime * r / u.^2 approximation
%          'gradient': gradient Richardson number [default]
%          'BV': Use Brunt-Vaisalla formulation with potential temperature
%
% OUT:
%  Ri0 = Source Richardson number at plume starting height (vent or water surface
%        level) for all runs

narginchk(1,2)
if nargin<2
    RiMode = 'simple';
end

if ~any([strcmp(RiMode,'simple')]) % strcmp(RiMode,'gradient') strcmp(RiMode,'BV')])
    RiMode = 'gradient';
    warning('Richardson: default to SIMPLE formula.')
end

g = 9.81;

nn = numel(qA.cI.Q);
Ri0 = nan(size(qA.cI.Q));

atmo = interpAtmoArray(qA.atmo,qA.cI.vh0(:));

% for ni = 1:nn

switch RiMode
    % Simple formulation, gprime * r / u.^2, typical for source experiments, but not great when the flow is not Boussinesq
    case 'simple' 
        rho_a0_stopgap = 1.12;
        gprime = g.* (rho_a0_stopgap - qA.pI.rho_B0)./max(cat(3,qA.pI.rho_B0, rho_a0_stopgap.*ones(size(qA.cI.Q))),[],3);
        Ri0 = gprime .* qA.pI.r_0 ./ qA.pI.u_0.^2;
  
  % --> GRADIENT AND BV MODES REQUIRE VERTICAL PROFILES TO GET GRADIENTS     
%     % More general gradient Richardson number
%     case 'gradient' 
%         Ri = g./ qA.pI.rho_B0 .* gradient(dat.pO.rho_B,dat.pO.z) ./ gradient(dat.pO.u,dat.pO.z);
% 
%     % Brunt-Vaisalla formulation (useful when using potential temperature 
%     % to get atmospheric conditions, but does not account for plume local density)
%     case 'BV'
%         atmo = interpAtmoArray(dat.pI.atmo,dat.pO.z);
%         N2 = getBruntVaisala(dat.pO.z,atmo(:,1),atmo(:,3),...
%     atmo(:,5)/100);
%         Ri = N2 ./ gradient(dat.pO.u,dat.pO.z);
        
end

% end


end