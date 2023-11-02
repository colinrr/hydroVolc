function phi = h2oExcessMass2ReservoirPorosity(gas_frac,cI)
% phi = h2oExcessMass2ReservoirPorosity(gas_frac,cI)
% Handy function to estimate phi0 for given Conduit starting conditions
    % gas_frac = mass fraction EXCESS gas above solubility that is present 
    %           (ie Csol * gas_frac = bulk_wt_%)
    % cI       = full getConduitSource structure (can be array)
    
    % Get vectors
    T   = [cI.T];
    dP  = [cI.dP];
    z0  = [cI.Z0];
    zW  = [cI.Zw];
    rhoMelt = [cI.rho_melt];
    
    rho_rock = 2400; % This should really be in getConSource
    rho_water = 1000;
    
    Pm = 9.81 .* ( rho_rock .* z0 + rho_water .* zW) + dP; % Reservoir pressure 
    Csol = meltH2Osolubility(Pm,0,T);
    
    n_h2o = Csol .* gas_frac; % Weight percent exsolved water
    
    % --- Add surface tension to get Pg? - test with vectors ---
    % '-> Recursive issue where we don't know r because we don't know phi - could solve, BUT...
    % '-> Tests give (2*ST/r) / Pm ~ 1e-3, so Pg~Pm is fine
    
%     load('fw_interpolant.mat','fw_interpolant_850C')
%     P.pm = Pm;
%     P.psat = Pm;
%     P.pb = findPb(P.pm,P.psat,T,fw_interpolant_850C);
%     ST = SurfaceTension(P,T);
%     
%     pg = Pm + 2*ST/r;
%     rhog = EoS_H2O_2(pg,Par.T);
%     mg = rhog * 4/3*pi*M3;
    % -------------
    
    [rho_h2o,~] = EoS_H2O_2(Pm,T); % Volatiles density - maybe also need surface tension to get pg??
    
    rhoBulk = (n_h2o./rho_h2o + (1-n_h2o)./rhoMelt).^(-1);
    
    phi = (rhoBulk - rhoMelt) ./ (rho_h2o - rhoMelt);
    
    assert(all(phi>=0),'Obtained negative porosity.')
    
end

function pb = findPb(pm,psat,T,fw_interpolant)
% Pressure in a bubble nucleus (Cluzel et al 2008)

KB = 1.38e-23;                      %Boltzman constant [J/K]

pb = zeros(size(pm));


for i = 1:length(pm)
    Vw = Molecular_VH2O(pm(i),T);
const = exp(Vw/KB/T*(pm(i)-psat(i)))...
    .* fw_interpolant(psat(i));

p1 = pm(i);
p2 = psat(i);

fun =  @(x) abs(fw_interpolant(x) - const); 
pb(i) = fminbnd(fun,p1,p2,optimset('TolX',1e-1)); 


end



end

function ST = SurfaceTension(Par,T)
% Hajimirza et al. JGR 2019

delta = .3200e-3;   % micro meter
alpha = .51;

Psat = Par.psat/1e6;    % Pa to MPa
Pb = Par.pb/1e6;
Pm = Par.pm/1e6;    % Pa to MPa
T = T - 273.15;         % K to C

ST_B = 1.21e-1 * exp(-2.24e-2*Psat)...
     + 1.47e-1 * exp(-1.90e-3*Psat)...
     + 7.5e-5 * (T-1000);
 
ST_inf = ST_B * (1-alpha);

ST = ST_inf + delta * (Pb-Pm);

end

function Vw = Molecular_VH2O(P,T)
% Molar volume of water (Ochs & Lange, Science 1999)

AV = 6.022e23;      % Avogadro number
T = T - 273.15;     % K to C
P =  P/1e5;         % Pa to bar


V = 22.9 + 9.5e-3*(T-1000) - 3.2e-4*(P-1);  % molar volume of water (cm^3)

Vw = V * 1e-6 / AV;

end
