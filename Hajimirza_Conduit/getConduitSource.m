function [conSource,atmo] = getConduitSource(varargin)
% conSource = getConduitSource(varargin)
% Get input struct for a run of Hajimirza/Gonnermann conduit model
% Inputs can be entered as a single struct or as Name/Value Pairs.
%
% CRowell Mar 2021

%% 1D Conduit DEFAULT input parameters

% ----- SOURCE PARAMS -----
Z0              = 6000;         % Initial depth (m)
dP              = 0;         % Initial overpressure (Pa)
T               = 850+273.15;          % Temperature (K)
theta           = 145;          % Contact angle between bubbles and crystals
ST_coeff        = 1/4*(2-cos(theta*pi/180))*(1+cos(theta*pi/180))^2;    % Heterogeneous nucleation factor
Q               = 1e8;          % Mass discharge rate (kg/s)
phi_frag        = .75;          % Critical porosity for fragmentation
N0              = 0;            % Initial bubble number density
phi0            = 0;            % Volume fraction of initial exsolved volatiles


conduit_radius  = 45;           % Conduit radius

composition     = struct();     % Struct containing optional name value pairs of chemical 
                                % WEIGHT FRACTIONS for Hui & Zhang 2007 viscosity model
                                % '-> see getComposition for default values
                                % field options:
                                %   'SiO2','TiO2','Al2O3','FeO','MnO','MgO','CaO','Na2O','K2O'
                                %
                                % NOTE! This provides some flexibility in composition, 
                                % but the water solubility model is still
                                % based on high-silica rhyolite, so use
                                % great caution (or get new
                                % solubility/diffusivity models) for SiO2
                                % much less than about 72-75 wt%
                                    
% Stuff exported from inside the condiut model
rho_melt    = 2400;       % Melt density [Kg/m^3]
rho_rock    = 2400;       % Average overburden rock density
% pf          = 1e5;        % Surface pressure - defaults to pressure from
                            % atmosphere file (based on vh0), but can be
                            % specified instead
                            
f0          = 0.0025;     % Friction coefficient

n0_excess       = 0;   % EXPERIMENTAL:
                        %  -> give a mass fraction of excess exsolved gas 
                        %     (relative to total DISSOLVED gas mass AT SATURATION)
                        %  -> ** REQUIRES non-zero N0 and phi0:
                        %      - N0 will be estimated internally if it is 0
                        %      - phi0 will be OVERWRITTEN since it is calculated
                        %      from n0_excess.


% Plume pressure coupling
vh0     = 00;              % vent altitude [m a.s.l.]
atmo = 'atmprofile.mat'; % default atmospheric profile. Input file must have [m x 5] array named 'atmprofile'
                         % Column order = [altitude temperature rel.hum. wind. pressure]
                         
Zw      = 0;                % water depth over vent [m]

% Potential exports to consider for later
% Par.K_melt      = 224e8;      % Bulk Modulus of melt
% Par.C0                        % Initial water content (calc'd from saturation)
% Par.mu          =             % Viscosity model could be adjusted to allow input for andesite/basalt/etc?
% Par.P_frag      = 1e5;        % What's that about? Apparently unused?

rho_l   = 1000; % Density of liquid water (kg/m^3)

% Conduit run failure thresholds to check for valid solution
zFailTol = 2;    % Z (depth) threshold in units of conduit radii: Zf/a < ZfailTol
mFailTol = 0.05; % Mach number tolerance for choking: abs(1-M) < MfailTol
pFailTol = .005; % Over-/underpressure threshold: (1 - Pm/Pf) < PfailTol

% Not conditional for valid solutions, but provides outcome interpretation
flareTol = .005; % Vent radius ratio tolerance to be considered flaring: abs(a/a_0 - 1) < flareTol

%% Parse input

    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    
    addParameter(p,'conduit_radius',conduit_radius)
    addParameter(p,'Z0',Z0)
    addParameter(p,'dP',dP)
    addParameter(p,'T',T)
    addParameter(p,'theta',theta)
    addParameter(p,'ST_coeff',ST_coeff)
    addParameter(p,'Q',Q)
    addParameter(p,'phi_frag',phi_frag)
    addParameter(p,'N0',N0)
    addParameter(p,'phi0',phi0)
    addParameter(p,'composition',composition)
    addParameter(p,'n0_excess',n0_excess)
    
    % Stuff ported out of original conduit script to allow input. Must be
    % added to Par struct internally
    addParameter(p,'rho_melt',rho_melt)
    addParameter(p,'rho_rock',rho_rock)
    addParameter(p,'f0',f0)
    addParameter(p,'pf',[])
    
    addParameter(p,'atmo',atmo)
    addParameter(p,'vh0',vh0)
    addParameter(p,'Zw',Zw)
   
    % Params for checking valid conduit solution after model run
    addParameter(p,'proxy',false)
    addParameter(p,'zFailTol',zFailTol)
    addParameter(p,'mFailTol',mFailTol)
    addParameter(p,'pFailTol',pFailTol)
    addParameter(p,'flareTol',flareTol)
    
    parse(p,varargin{:})
    conSource = p.Results;
    clear atmo;
    
    % Populate chemical composition
    conSource.composition = getComposition(conSource.composition);
    
    % Atmospheric profile
    if(ischar(conSource.atmo))
        load(conSource.atmo)
        if or(istable(atmprofile),isstruct(atmprofile))
            [atmo.atmo,atmo.atmoVars,atmo.atmoUnits] = getAtmoArray(atmprofile);
        else
            atmo.atmo = atmprofile(~isnan(atmprofile(:,1)),:);
            atmo.atmoVars = {'Altitude', 'Temperature', 'Pressure', 'windAbs', 'relativeHumidity'};
            atmo.atmoUnits = {'m a.s.l.', 'K', 'Pa','m/s', '%'};
        end
%     elseif(isnumeric(plumeSource.atmo))
%         assert(size(plumeSourc.atmo,2)==5,'Atmo profile must be [m x 5]')
    end
    assert(size(atmo.atmo,2)==5,'Atmo profile must be [m x 5]')

    % Get ambient pressure at vent - default from atmo profile
    if isempty(conSource.pf)
        conSource.pf = 10.^interp1(atmo.atmo(:,2),log10(atmo.atmo(:,1)),conSource.vh0,'pchip','extrap') + rho_l*9.81*conSource.Zw;
    end
    
%     conSource = struct2table(conSource,'AsArray',true); % table can be deprecated

    puf = fieldnames(p.Unmatched);
    if ~isempty(puf)
        warning('Unrecognized 1Dconduit Name/Value argument(s):')
        fprintf('\t%s\n',puf{:})
    end
    
    % Run some checks
    if conSource.n0_excess ~= 0
        if ~ismember('phi0',p.UsingDefaults)
            warning('Input var "phi0" will be overwritten by n0_excess calculation.')
        end
        conSource.phi0 = h2oExcessMass2ReservoirPorosity(conSource.n0_excess,conSource);
        
        % Estimate N0 if needed - note that this method is very approximate
        % and only tested for a very limited rhyolitic composition
        if ismember('N0',p.UsingDefaults)
            conSource.N0 = linearN0estimator(conSource.Q,conSource.pf);
        end
    end
    
    if conSource.phi0==0 && conSource.N0~=0
        warning('Cannnot have non-zero Bubble Number Density with zero initial porosity. BND set to 0.')
        conSource.N0 = 0;
    end
    
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