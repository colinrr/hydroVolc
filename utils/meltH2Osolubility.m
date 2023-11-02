function [Cw,Cc] = meltH2Osolubility(P,Xc,T)
% Solubility of water and CO2 in melt
% Liu et al 2005
%   P = Pressure [Pa]
%   Xc = mole fraction of CO2 in total volatile content (assuming H2O+CO2)
%   T  = Temperature
%
% Exported rom Hajimirza conduit model.

    Xw = 1 - Xc;
    Pw = Xw .* P/1e6;
    Pc = Xc .* P/1e6;
  
    a1 = 354.94;
    a2 = 9.623;
    a3 = -1.5223;
    a4 = 0.0012439;
    a5 = -1.084e-4;
    a6 = -1.362e-5;
    %
    b1 = 5668;
    b2 = 0.4133;
    b3 = 2.041e-3;
    b4 = -55.99;
    
    
    Cw = (a1*Pw.^(1/2) + a2*Pw + a3*Pw.^(3/2))./T + a4*Pw.^(3/2) +...
      Pc.*(a5*Pw.^(1/2) + a6*Pw); %H2O content in wt.%
    Cc = b1*Pc/T + Pc.*(b2*Pw.^(1/2) + b3*Pw.^(3/2)) + b4*Pc.*Pw/T;
  
    Cw = Cw/100; % fraction
    Cc = Cc/1e6; % fraction
    
end