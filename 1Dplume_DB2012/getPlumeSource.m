function pI = getPlumeSource(varargin)
% plumeInput = getPlumeSource(varargin)
% Get input struct for a run of Aubry/Jellinek 2018 1D plume model 
% (modified from Degruyter and Bonadonna 2012).
% Inputs can be entered as a single struct or as Name/Value Pairs.
%
% C Rowell Mar 2021

%% 1D plume DEFAULT input parameters

% ----- SOURCE PARAMS ----- (original)
T0      = 1150;     % exit temperature                  [K]
n_0     = 0.05;     % TOTAL water mass fraction         [wt.%]
vh0     = 00;       % vent altitude                     [m a.s.l]
u_0     = 200;      % exit velocity                     [m/s]
r_0     = 200;      % vent radius                       [m]
xv_0    = 1;        % steam/water dryness fraction      -

% ----- MAGMA PROPS ----- (exported from hmodel.m)
rho_m   = 2400;     % magma density                     [kg/m^3]
C_s     = 1250;     % specific heat of solid pyroclasts [J/kg/K]

% ----- Default GSD params -----
phi0        = 0.75; % Large pyroclast porosity
D           = 3.0;  % GSD Power law exponent
phiSz_min   = -9;   % GSD min phi size (max r)
phiSz_max   = 10;   % GSD max phi size (min r)
Pg          = 1e6;   % Bubble gas pressure. Should really inherit from conduit, or perhaps assume atmospheric

% The next 5 params define the GSD and will be generated from the above 5
% by default. See getGSD for more info.
nsi         = [];   % particle mass fraction at each grain size (including bubble gas). Sum(mi) = 1.
ni          = [];   % mass fraction of bubbles at each grain size, relative to particle mass fraction.
Rgsd        = [];   % particle radii (m)
rhoi        = [];   % bulk pyroclast density as function of grain size (kg/m^3)
pori        = [];   % pyroclast porosity as a function of grain size

% ----- COEFFICIENTS ------
lambda  = 10^(-2);  % lambda    condensation rate           [s^-1]
alpha   = 0.1;      % alpha     radial entrainment coeff.   -
beta    = 0.5;      % beta      wind entrainment coeff.     -

% ----- ATMO PROFILE -----
atmo      = 'atmprofile.mat'; %   [m x 5] table or array, m = height levels, with columns:
atmoVars  = {'Altitude', 'Temperature', 'Pressure', 'Wind_abs', 'Relativehumidity'};
atmoUnits = {'m a.s.l.', 'K', 'Pa','m/s', '%'};

    % NOTE: Tables may alternatively have variables: 
    %       'Meridionalwindspeed', 'Zonalwindspeed'  - [m/s]

% ----- PLUME ENTRAINMENT MODEL SWITCHES -----
wind = true; % t/f to turn wind effects on/off. Default is true (include wind).
% baratio   beta/alpha ratio for entrainment model 4
model = 2;  % entrainment model to use (see Aubry and Jellinek EPSL 2018)
nexp  = 1;  % normalization exponent as in Devenish et al (see Aubry and Jellinek EPSL 2018)


% ----- OTHER ------
% Parameters for model coupling, decompression, MWI, etc that are not
% directly needed for the plume model
% rho_B0  	source bulk density
% rho_g0    source gas density
% Ld      	decompression length

% ----- MAGMA-WATER INTERACTION PARAMS -----
mwiModel = 2; % MWI model version: 1 = equilibrium mixing model w/ static PSD (out-dated)
                                 % 2 = equilibrium mixing model w/ dynamic PSD
                                 % 3 = pyroclast-water disequilibrium model (in-dev)
Es           = 100; % surface energy density of pyroclasts   [J/m^2]
Tw0          = 274.15; % Surface water temperature           [K]
saScale      = 10;  % Particle surface area roughness scale (multiplied w/ spherical surface area)
fragPar      = 0.1; % Fragmentation energy partitioning fraction.
alphaRTscale = 0; % [0-1] Scale factor f for Rayleigh-Taylor water entrainment. 
                    % '-> Experimental. alpha = (1-f)*alpha_Ri + f*alpha_RT

k_m     = 1.5;      % melt thermal conductivity         [W/m/K] (MWIv3 only)
T_g     = 700;      % Glass transition temperature minimum [K] (normally set inside MWI model, enter a value to override)
T_g_rng = 50;      % Range of glass transition temperatures [K]

% MWI ENTRAINMENT SCALES
useJetEntranceLength = true;  % Scale entrainment coefficient to a jet entrance length
useDecompressLength  = true;  % Prevent entrainment over a distance defined by vent overpressure

% For MWIv2: fragmentation output grain size mean and standard dev. (phi units)
% Default numbers correspond to Askja Unit C from Costa et al. (2016) "Assessing total grain size..."
phiFrag_mu   = 3.43; 
phiFrag_sig  = 1.46;
phiFrag_cutoff = [];
k_w     = 0.6; % Water thermal conductivity

% --------------------------------------------

%% Parse input
%     narginchk(0,1)

    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    
    addParameter(p,'T0',T0)
    addParameter(p,'n_0',n_0)
    addParameter(p,'vh0',vh0)
    addParameter(p,'u_0',u_0)
    addParameter(p,'r_0',r_0)
    addParameter(p,'xv_0',xv_0)
    
    addParameter(p,'lambda',lambda)
    addParameter(p,'alpha',alpha)
    addParameter(p,'beta',beta)

    addParameter(p,'atmo',atmo)
    addParameter(p,'atmoVars',atmoVars)
    addParameter(p,'atmoUnits',atmoUnits)

    addParameter(p,'wind',wind)
    addParameter(p,'baratio',NaN)
    addParameter(p,'model',model)
    addParameter(p,'nexp',nexp)
    
    addParameter(p,'rho_m',rho_m)
    addParameter(p,'C_s',C_s)
    addParameter(p,'k_m',k_m)
    addParameter(p,'T_g',T_g)
    addParameter(p,'T_g_rng',T_g_rng)
    
    % GSD
    addParameter(p,'phi0',phi0)
    addParameter(p,'phiSz_min',phiSz_min)
    addParameter(p,'phiSz_max',phiSz_max)
    addParameter(p,'D',D)
    addParameter(p,'Pg',Pg)
    addParameter(p,'nsi',nsi)
    addParameter(p,'ni',ni)
    addParameter(p,'Rgsd',Rgsd)
    addParameter(p,'rhoi',rhoi)
    addParameter(p,'pori',pori)
    
    % MWI
    addParameter(p,'mwiModel',mwiModel)    
    addParameter(p,'Es',Es)
    addParameter(p,'Tw0',Tw0)
    addParameter(p,'k_w',k_w)
    addParameter(p,'saScale',saScale)
    addParameter(p,'fragPar',fragPar)
    addParameter(p,'alphaRTscale',alphaRTscale)    
    addParameter(p,'useJetEntranceLength',useJetEntranceLength)
    addParameter(p,'useDecompressLength',useDecompressLength)
    addParameter(p,'phiFrag_mu',phiFrag_mu)
    addParameter(p,'phiFrag_sig',phiFrag_sig)
    addParameter(p,'phiFrag_cutoff',phiFrag_cutoff)
    
    % Other
    addParameter(p,'rho_g0',[])
    addParameter(p,'rho_B0',[])
    addParameter(p,'Ld',[])
%     addParameter(p,'wScale',[])
%     addParameter(p,'m_ew0',0)
%     addParameter(p,'mwe',NaN)
%     addParameter(p,'dD',0)

    parse(p,varargin{:})
    
    pI = p.Results; % Plume Source struct
    
    % --- Atmospheric profile ----
%     if isempty(plumeSource.atmo)
%         load('atmprofile.mat')
%         plumeSource.atmo = atmprofile;
    if(ischar(pI.atmo))
        load(pI.atmo)
        if or(istable(atmprofile),isstruct(atmprofile))
            [pI.atmo,pI.atmoVars,pI.atmoUnits] = getAtmoArray(atmprofile);
        else
            pI.atmo      = atmprofile(~isnan(atmprofile(:,1)),:);        
        end
%     elseif(isnumeric(plumeSource.atmo))
%         assert(size(plumeSourc.atmo,2)==5,'Atmo profile must be [m x 5]')
    end
    assert(size(pI.atmo,2)==5,'Atmo profile must be [m x 5]')
    
    % ---- Generate GSD if needed ----
    if isempty(pI.nsi)
        X = 1;
        [pI.nsi,pI.ni,pI.Rgsd,pI.rhoi,pI.pori] = getPSD(pI.rho_m, pI.phi0, pI.Pg, pI.T0, 1, pI.D, pI.phiSz_min, pI.phiSz_max);
        
        % Re-calculate free gas fraction    
        % NOTE! Assumes n_0 is total exsolved gas and incorporates some into porosity
        pI.n_0 = pI.n_0 - (1-pI.n_0).*sum(pI.ni.*pI.nsi);
    end
    
    % Unmatched fields?
    puf = fieldnames(p.Unmatched);
%     if ~isempty(puf)
%         warning('Unrecognized 1Dplume Name/Value argument(s):')
%         fprintf('\t%s\n',puf{:})
%     end
    
end

% function [A,vars,units] = getAtmoArray(atmtable)
%     
%     vars = atmtable.Properties.VariableNames;
%     atmtable = rmmissing(atmtable);
%     units = atmtable.Properties.VariableUnits;
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