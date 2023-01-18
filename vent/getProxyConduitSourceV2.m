function [cI,cO,atmo] = getProxyConduitSourceV2(varargin)

% Get structs for vent parameters for input into MWI model, without running
% the conduit model.
%    Specialized optional inputs for external water infiltrating the conduit:
%    (ie unique to this function, not present in conduit model)
%       n_ec = mass fraction of external water
%       T_ec = temperature of external water
%       Z_ec = average depth of aquifer
%       Cp_s = Magma heat capacity
%       Q0   = magma discharge rate (not counting external water)
%
%    Calculated interally and used to approximate T after water mixing:
%       Cp_v = water vapor heat capacity
%       L_v  = Latent heat of vaporization
%       Cp_w = liquid water heat capacity
%       Tr   = reference temperature for magma heat capacity
%
% Inputs can be entered as a single struct or as Name/Value Pairs.
%
% CRowell
% - V2.0, Jan 2022. No choking?. ??

%% Conduit Proxy DEFAULT input parameters

T           = 850+273.15; % Temperature (K)
rho_melt    = 2400;       % Melt density [Kg/m^3]
K_melt      = 224e8;      % Bulk Modulus of melt (Pa)
Cm          = 0.005;       % Residual dissolved water mass fraction
n_0         = .037;       % Free gas mass fraction
a           = 45;         % Conduit radius (m)
g           = 9.81;       % m/s^2
Cp_s        = 1250;       % Pyroclast heat capacity (J/kg/K)
rho_e       = 1000;       % Surface water density
R_v         = 461;        % gas constant of volcanic gas (water) (J/kg/K) 
phi_frag    = .75;        % Critical porosity for fragmentation

% Plume pressure coupling
vh0     = 00;              % vent altitude [m a.s.l.]
atmo    = 'atmprofile.mat'; % default atmospheric profile. Input file must have [m x 5] array named 'atmprofile'
                         % Column order = [altitude temperature rel.hum. wind. pressure]
                         
Zw        = 0;             % water depth over vent [m]
beta_crit = 0.95;         % Critical overpressure value, below which vent is assumed not choked

Q        = 1e8;           % Initial magma mass discharge rate (kg/s)

% Conduit water params
n_ec     = 0;             % mass fraction of external water in the conduit: q_ec/(Q0+q_ec)
rho_rock = 2400;          % Rock density
T_ec     = 275;           % External water temp (K)
Z_ec     = 100;           % Average aquifer depth (m) (for water props)
% Tr       = 274.15;        % Reference temperature

%% Parse input

    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    
    addParameter(p,'conduit_radius',a)
    addParameter(p,'T',T)
    addParameter(p,'rho_melt',rho_melt)
    addParameter(p,'K_melt',K_melt)
    addParameter(p,'Cm',Cm)
    addParameter(p,'Cp_s',Cp_s)
    addParameter(p,'n_0',n_0)
    addParameter(p,'Q',Q)
    addParameter(p,'atmo',atmo)
    addParameter(p,'vh0',vh0)
    addParameter(p,'Zw',Zw)
    addParameter(p,'pf',[])
    addParameter(p,'n_ec',n_ec)
    addParameter(p,'T_ec',T_ec)
    addParameter(p,'Z_ec',Z_ec)
    addParameter(p,'rho_rock',rho_rock)
    addParameter(p,'phi_frag',phi_frag)
    addParameter(p,'beta_crit',beta_crit)
    addParameter(p,'conMod','proxy') % 
%     addParameter(p,'Tr',Tr)

    parse(p,varargin{:})
    cS = p.Results;
    
    clear atmo;
%% Calculate params

    % ----- Atmospheric profile ------
    if(ischar(cS.atmo))
        load(cS.atmo)
        if or(istable(atmprofile),isstruct(atmprofile))
            [atmo.atmo,atmo.atmoVars,atmo.atmoUnits] = getAtmoArray(atmprofile);
        else
            atmo.atmo = atmprofile(~isnan(atmprofile(:,1)),:);
            atmo.atmoVars = {'Altitude', 'Temperature', 'Pressure', 'Wind_abs', 'Relativehumidity'};
            atmo.atmoUnits = {'m a.s.l.', 'K', 'Pa','m/s', '%'};
        end
%     elseif(isnumeric(plumeSource.atmo))
%         assert(size(plumeSourc.atmo,2)==5,'Atmo profile must be [m x 5]')
    end
    assert(size(atmo.atmo,2)==5,'Atmo profile must be [m x 5]')

    % Get ambient pressure at vent - default from atmo profile
    if isempty(cS.pf)
        cS.pf = 10.^interp1(atmo.atmo(:,2),log10(atmo.atmo(:,1)),cS.vh0,'pchip','extrap') + rho_e*g*cS.Zw;
    end
    

    % ----- CONDUIT EXTERNAL WATER and VENT PARAMS ---------
    A_c     = pi * cS.conduit_radius.^2;
%     cS.T    = cS.T+273.15;
    cS.Q0   = cS.Q;
    
    if cS.n_ec~=0
        cS.Q    = cS.Q./(1-cS.n_ec);
        
        % Adjusted gas mass fractions
        n0_pr   = cS.n_0.*(1-cS.n_ec);
        n       = cS.n_ec + n0_pr;
        
        % Water props
        Pw      = g*cS.Z_ec*cS.rho_rock + cS.pf  ;
        Wec      = waterProps(Pw,cS.T_ec);
        
        % Rough estimates for water latent heat and heat capacities
        L_ec           = Wec.h_vs - Wec.h_ws; % Latent heat estimate
        Tv = mean([Wec.Tsat cS.T]);
        Tw = mean([cS.T_ec Wec.Tsat]);
        rhow = density(Pw,Tw);
        rhov = density(Pw,Tv);
        Cp_w = heatCapacityP(rhow,Tw);
        Cp_v = heatCapacityP(rhov,Tv);
        Cb_pr = (1-n).*cS.Cp_s + n0_pr*Cp_v; % Weighted heat cap. of initial mixture
        
        % Mixture temparature estimate
        Tf    = (cS.T.*Cb_pr - cS.n_ec.*Cp_w.*(Wec.Tsat - cS.T_ec) - cS.n_ec.*L_ec + cS.n_ec.*Cp_v.*Wec.Tsat) ./ (Cb_pr + cS.n_ec.*Cp_v);

        % Vent pressure estimate w/ 1/.86 empirical factor due to differing
        % sound speed calc in the model
%         P_v     = (n.*R_v.*Tf).^(1/2).*cS.Q./A_c;
        P_v     = 1/0.97*(n.*R_v.*Tf).^(1/2).*cS.Q./A_c;
        
        
    else
%         cS.Q    = cS.Q0;
        n       = cS.n_0;
        Tf      = cS.T;
%         P_v      = cS.Q./A_c.*(n.*R_v.*cS.T).^(1/2);
        % 0.97 factor scales to match model thresholds for Mach #
        P_v      = 1/0.97*cS.Q./A_c.*(n.*R_v.*cS.T).^(1/2);
    end
    
    % Check vent overpressure
    if P_v/cS.pf <= cS.beta_crit
        P_v = cS.pf;
        cO.Par.choke = 0;
    else
        cO.Par.choke = 1;
    end
    
    % Recalculate vent parameters with estimated pressure
%     W = waterProps(P_v,Tf);
    [cO.rho_g,K_g] = EoS_H2O_2(P_v,Tf);  % Gas density and bulk mod
    cO.rho_magma   = (n./cO.rho_g + (1-n)./cS.rho_melt).^(-1); % Bulk density
    cO.porosity    = n.*cS.rho_melt ./ (n.*cS.rho_melt + (1-n).*cO.rho_g); % Total porosity include bubble fraction
    K_magma = 1./((cO.porosity./K_g) +...     % Mixture bulk modulus
        (1-cO.porosity) ./ cS.K_melt);
    cO.C = (K_magma ./ cO.rho_magma).^(1/2);  % Mixture sound speed
    
%% Populate output structs
 
% 1D Conduit input parameters
cI.T        = cS.T; % Tf? - check to make sure this isn't used later
cI.atmo     = cS.atmo;
cI.vh0      = cS.vh0;
cI.Zw       = cS.Zw;
cI.conduit_radius = cS.conduit_radius;
cI.n_0      = cS.n_0;
cI.n_ec     = cS.n_ec;
cI.pf       = cS.pf;
cI.phi_frag = cS.phi_frag;
cI.Q        = cS.Q0;

% 1D Conduit output parameters
cO.Par.proxy    = true;
cO.Par.T        = Tf;
cO.Par.rho_melt = cS.rho_melt;
cO.Par.Q        = cS.Q;   % New total mass flux after incorporation of water
cO.Par.Q0       = cS.Q0;  % Initial magma+volatile mass flux
cO.Par.frag     = 1;
cO.Par.pf       = cS.pf;
cO.Par.K_melt   = cS.K_melt;
cO.Par.phi_frag = cS.phi_frag;
cO.Par.Zf       = NaN;
cO.Par.conValid = NaN;
cO.Par.c_s      = cS.Cp_s;
cO.Par.rho_rock = cS.rho_rock;
cO.Par.T_ec     = cS.T_ec;
cP.Par.Z_ec     = cS.Z_ec;
cO.pg           = P_v;
cO.M3           = cS.phi_frag/(1-cS.phi_frag)*3/(4*pi);
cO.Cm           = cS.Cm;
% √ cO.porosity
% √ cO.rho_g        
% √ co.rho_magma
cO.U = cO.Par.Q./(cO.rho_magma.*A_c);
cO.n = n;
cO.pm = cO.pg;
cO.K = cO.pm./cI.pf;
cO.Z = 0;
cO.M = cO.U/cO.C;
cO.cin = NaN;

% √ cO.C
cO.a = cI.conduit_radius;

cI.proxy = true;

end

% function [A,vars,units] = getAtmoArray(atmtable)
%     
%     vars = atmtable.Properties.VariableNames;
%     units = atmtable.Properties.VariableUnits;
%     atmtable = rmmissing(atmtable);
%     A = zeros(size(atmtable,1),5);
%     A(:,1) = atmtable.Pressure;
%     A(:,2) = atmtable.Altitude;
%     A(:,3) = atmtable.Temperature;
%     if ismember('Wind_abs',vars)
%         A(:,4) = atmtable.Wind_abs;
%     else
%         wU     = atmtable.Meridionalwindspeed;
%         wV     = atmtable.Zonalwindspeed;
%         A(:,4) = (wU.^2 + wV.^2).^(1/2); 
%     end
%     A(:,5) = atmtable.Relativehumidity;
%     
% end