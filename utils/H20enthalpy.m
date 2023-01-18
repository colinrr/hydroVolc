function [Hv,Hl,Hi,Hd] = H20enthalpy(T)
 % Returns enthalpies of water phases as a function of
 % temperature, as well as dry air.
 %
 % Based on Mastin 2007
 % C Rowell May 2021
 
 assert(isvector(T),'T must be a vector')
 if size(T,2)>1
     T = T';  % Needs to be a column vector
 end
 
    R_d = 287;      % gas constant of dry air (J/kg/K) 
    R_v = 461;      % gas constant of volcanic gas (water) (J/kg/K)
    Hv0 = 2.5007e6; % Reference vapour enthalpy at 273.15
    Hd0 = 2.73e6;   % Reference dry air enthalpy at 270 K
    Hi0 = -333430;  % Reference ice enthalpy at 273.15 K
    Ci  = 1850;     % Ice heat capacity

    Ti2 = 273.15-40;
    Ti1 = 273.15-17;
    Td0 = 270;
    Tf  = 273.15;    % Freezing temp [K]
    Tb  = 373.15;    % Boiling temp [K]
    Tc  = 647.25;    % Critical temp [K]

    % ai = [-7.8889166 2.5514255 -6.716169 33.239495 105.38479 174.35319 148.39348 48.631602]';
    % Enthaply coefficients for vapour and dry air
    av = [4.070 -1.108e-3 4.152e-6 -2.964e-9 8.07e-11];
    ad = [3.653 -1.337e-3 3.294e-6 -1.913e-9 2.763e-11];

    hl_vec = [% [K]    [J/kg]  
                273.25 381.14
                283.25 42406
                293.25 84254
                303.25 126090
                313.25 167920
                323.25 209750
                333.25 251570
                343.25 293430
                353.25 335350
                363.25 377350
                373.25 419490];

   
    % Need to create a T switch here?
%     switch T

    % Heat capacities as f(T)
    CvFun = @(T) R_v * sum(av'.*T.^((0:4)'));
    CdFun = @(T) R_d * sum(ad'.*T.^((0:4)'));
    
    Hv = Hv0 + integral(CvFun,Tf,T);
    Hd = Hd0 + integral(CdFun,Td0,T);
    Hl = interp1(hl_vec(:,1),hl_vec(:,2),T,'linear','extrap');
    Hi = Hi0 + Ci*(T-Tf);

end