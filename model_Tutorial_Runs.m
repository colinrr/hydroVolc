
%% Test script of coupled conduit/plume model, single run
%
% CRowell Mar 2021
%
clear all; close all

% Run from this script's own directory if you use this vv
addpath(genpath('.')) % Consider using/defining absolute paths for data etc

%% INPUT EXAMPLES - comment/uncomment as you like, 
which_scenario = 1; % Run which case?

% I give two test cases for each, with slightly different input. Note that
% I've set water depths to values that keep things running quickly. As you
% approach Zw ~ (5 x conduit_radius), things will slow down a lot.
switch which_scenario
    case 1
%  ###################### (1) PROXY RUN, KATLA #######################
    % ---- CONDIUT  ----
    cI      = [];
    cI.vh0  = 850;                % vent altitude (m)
    cI.Q    = 5e7;                  % magma mass discharge (kg/s)
    % vv Use an empirical fit function from running many conduit models to get vent
    % radius (handy to give an idea of radius values that are appropriate)
    cI.conduit_radius = extrapVentRadius(cI.Q); 
    cI.atmo = 'atmoFiles/Katla_avg_atmo_Oct12_1979-2005.mat';
    cI.rho_melt = 2600;   % Melt density
    cI.proxy    = true;   % Proxy conduit run?
    cI.T_ec     = 273.15; % external GROUNDwater temp
    cI.rho_rock = 2800;   % lithostatic density
    cI.T        = 1500;   % magma Temp (K) - [1470 to 1540]; % Basalt, 1200-1270 C
    cI.n_0      = 0.01;   % Reasonable gas (H2O) mass fraction for basalt
    cI.Zw       = 20;     % Surface water depth (m)
    cI.n_ec     = 0.0;    % Prescribed mass fraction of groundwater infiltration into the conduit
    % ------ MWI/PLUME ------
    pI.Tw0      = cI.T_ec;  % external SURFACE water temp
    pI.useDecompressLength = false; % Turn off explosive decompression entrainment limiter
                                % (note at present this introduces a somewhat
                                % non-physical expansion at the vent (just
                                % happens instantly, essentially)

    % Glass transition data for fragmentation:
    % Giordano et al 2005, basalt assuming ~0.5wt% residual H20
    pI.T_g       = 870; % K
    pI.T_g_rng   = 50;
    pI.phiFrag_mu   = 3; % Very rough eyeballed output sizes for second mode of Jonsdottir 2015
    pI.phiFrag_sig  = 1.25;
    pI.D         = 2.6; % PSD exponent - initially very coarse particle size distribution (very fine would be about 3 or 3.2)

    % Second run - duplicate, but change water depth
    cI2 = cI;
    cI2.Zw = 40;
    pI2 = pI;

case 2
%  ############ (2) HUNGA-TONGA BIG ERUPTION, 2 MER values #############
    cI.proxy    = true;
    cI.T_ec     = 302; % ERA5 SST
    cI.rho_rock = 2800;
    cI.atmo     = 'atmoFiles/HTHH_04_2022-01-15.mat';
    cI.T        = 1173; % Andesite
    cI.n_0      = 0.03; % .015 to 0.05?
    cI.n_ec     = 0.1;
    cI.Q        = 8e8; % She went big - 8e8 to 5e9? (number I actually put in below is probably pretty unreasonable)
    cI.conduit_radius = extrapVentRadius(cI.Q); 

    cI.Zw       = 150;

    % ------ MWI/PLUME ------
    pI.Tw0      = cI.T_ec;
    % pI.D        = 2.8; % A bit coarser for basaltic-andesite?
    % Glass transition data - Colombier et al 2018;
    pI.T_g       = 808;
    pI.T_g_rng   = 50;

    cI2 = cI;
    cI2.Q = 10e9; % Volcano make real big boom 
    cI2.conduit_radius = extrapVentRadius(cI2.Q); % Need to roughly adjust radius to match Q
    pI2 = pI;
    

case 3
%  ############ (3) CONDUIT RUN, RHYOLITE PLINIAN, 2 latitudes #############
    % Note that running the conduit model is a bit tricky - solutions are not
    % "automatically" physical, but must be found (ie a "shooting"
    % solution). That is, the fragmented flow must either reach the choking 
    % condition (Mach 1 and overpressured) almost exactly on the vent, OR 
    % be less than Mach 1 and pressure-balanced. 
    % See "conduitRadiusFromQ" and "conduitQfromRadius" for
    % functions that enable this. Adding chamber overpressure etc would
    % require similar.

    cI.conduit_radius = 45.45;
    cI.Q = 1e8;
    cI.atmo = 'atmoFiles/atm_ERAreanalysis_Grimsvotn2011_01_absWind.mat';
    cI.Zw = 1500;
    
    pI.saScale   = 50; % Very high particle surface area because why not
    pI.wind      = false; % Can play with wind
    
    cI2 = cI;
    cI2.dP = 1e6; % Add some magma chamber overpressure to mess with it a little
    pI2 = pI;
    cI2.atmo = 'atmoFiles/atm_ERAreanalysis_Tungarahua2014_01_absWind.mat';
    
end





%% DO THE THING 
tic
dat1 = runCoupledModel(cI,pI);
toc;tic
dat2 = runCoupledModel(cI2,pI2); %,[],20);
toc

% Summary table - note this function is out of date but deserves an (easy) update
% singleRunSummary(dat1)
% singleRunSummary(dat2)

%% Plot


% These versions work but are optimized for manuscripts. Start by changing
% font size, haha ('fs' in the functions). Also I should standardize the 
% input better, not all of these will work for every case...
if which_scenario==3
    plotConduitPair_papV2(dat1,dat2)
end

if ~dat1.dO.Ld_gt_Zw % have to check if the MWI model actually ran
    mwiax = plotMWIoutput_pap(dat1.wO,dat1.cI);
else
    mwiax = []; 
end
if ~dat2.dO.Ld_gt_Zw
    plotMWIoutput_pap(dat2.wO,dat2.cI,mwiax);
end

pause(0.5) % Plotting goes weird if you rush it sometimes...
if which_scenario==3
    [aax] = plotAtmoProfile(cI.atmo);
    plotAtmoProfile(cI2.atmo,aax);
end

plotPlumePairMS([dat1 dat2]); % not optimized for different atmospheres so watch out for the tropopause location