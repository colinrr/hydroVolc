% Running some tests - CR 2021

P = 2e6;
T = 273.15+850;

% Enthalpy
rho = density(P,T);    % kg/m^3
h   = enthalpy(rho,T); % J/kg

% Saturation enthalpy
Tsat = saturationTemperature(P);
rho_ws = density(P,Tsat-1e-5);
rho_vs = densityTX(Tsat,1);
[hw,~] = enthalpy(rho_ws,Tsat);
[hv,~] = enthalpy(rho_vs,Tsat);

% Wet steam mixture density at saturation T with dryness fraction x
X = 0.5;
rho_x = densityTX(Tsat,X);

%% 
m = 1e8;
n = 0.03;
u = 150;
rho = 30;
r = (m/(pi*u*rho)).^(1/2);
rho_a = 1000;
alph = 0.05;
T = 273.15+850;
Tw0 = 274.15;
P = 5e5;
F = 0.9;

dm_w_dz = 2 * pi * r * u * alph * sqrt(rho*rho_a);

W = waterProps(P,T);
W0 = waterProps(P,Tw0);
dEfrag = dm_w_dz * (W.h - W0.h) / (F * (1-n) * m)
