function [a, a_by_mass] = estimateConduitRadius(cI)
    % Crude calcs to throw out rough range estimates for a

    % rough bounds on estimated values
    % start pre-frag post-frag
    g = 9.81;
    u0 = [2 5 10];
    u   = [10 40 150];
    du  = [diff(u(1:2)) .* [1 1] diff(u(2:3))];
    rho = [cI.rho_melt 400 100];
    f0 = 0.0025;
    
    mu0 = viscosity(cI.T,meltH2Osolubility(cI.rho_rock * cI.Z0 * 9.81 + cI.dP + cI.pf,0,cI.T),cI.composition);
    mu  = [mu0 2.4e6 1e-3];
    
    % a by mass
    a_by_mass = (cI.Q./(rho(1) .* u0 .* pi)).^(1/2);
    
    
    DeltaP = cI.rho_rock * cI.Z0 * g + cI.dP - cI.pf;
    
    Delta = DeltaP + u .* du + rho .* g .* cI.Z0;
    
    a = [(8 .* cI.Z0 .* mu(1:2) .* u(1:2) ./ Delta(1:2)).^(1/2) ...
        f0 .* cI.Z0 .* rho(3) .* u(3).^2 ./ Delta(3)] ;
    
    Delta_mean = DeltaP + mean(u(1:2)) .* du(1) + mean(rho(1:2)) .* g .* cI.Z0;
    a_by_mean = (8 .* cI.Z0 .* mean(mu(1:2)) .* mean(u(1:2)) ./ Delta_mean).^(1/2);
    
end