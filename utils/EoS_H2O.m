function [Pg, Kg] = EoS_H2O(rho,T)
% Exported from Hajimirza conduit model for easy use - CR Mar 2021
% Modified Redlich and Kwong EoS for water vapor from Holloway 1977


R = 83.12;          % Gas constant cm^3.bar/(deg mole)
M = 18.01528e-3;    % Molar mass of water kg/mol

TC = T - 273.15;    % Degree Cel

ao = 35e6;
b = 14.6;
a = 166.8e6 - 193080*TC + 186.4*TC.^2 - 0.071288*TC.^3;



V = (1./rho)*M*1e6;      % molar volume (cm^3/mol)

Pg = R*T./(V-b) - a./(sqrt(T).*V.*(V+b));

Pg = Pg * 1e5;          % bar to Pa


% Bulk Modulus
dpdv = -R*T./(V-b).^2 + a./(sqrt(T)) * (2*V+b)./(V.*(V+b)).^2;
dvdrho = -(1./rho).^2 * M * 1e6;
Kg = 1.33 * rho .* dpdv .* dvdrho; % 1.33 is Cp/Cv
Kg = Kg * 1e5;

Kg(Pg==0) = inf;

end

