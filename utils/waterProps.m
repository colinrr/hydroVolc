function W = waterProps(P,T)
% Returns water and steam properties for pressure P (Pa) and temperature T
% (K). Uses IAPWS95 MATLAB code from Junglas, P. 2009

    % Properties at P,T
    W.rho = density(P,T);    % kg/m^3
    W.h   = enthalpy(W.rho,T); % J/kg
    W.C_p       = heatCapacityP(W.rho,T);
    W.C_v       = heatCapacityV(W.rho,T);
    W.gamma     = W.C_p/W.C_v;                  % Ratio of specific heats
    W.c         = soundVelocity(W.rho,T);
    W.Kg        = W.c.^2 * W.rho;               % Adiabatic bulk modulus

    % Properties at saturation temp for steam and water
    W.Tsat      = saturationTemperature(P);
    W.rho_ws    = density(P,W.Tsat-1e-5);
    W.rho_vs    = densityTX(W.Tsat,1);
    W.h_ws      = enthalpy(W.rho_ws,W.Tsat);
    W.h_vs      = enthalpy(W.rho_vs,W.Tsat);
    
    
end