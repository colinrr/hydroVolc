
%% Test script of coupled conduit/plume model, single run
%
% CRowell Mar 2021
%
clear all; close all

% addpath(genpath('./Hajimirza_Conduit/'))
% addpath(genpath('./1Dplume_DB2012/'))
% addpath(genpath('./vent/'))
% addpath(genpath('./utils/'))
% addpath(genpath('./plot-tools/'))
addpath(genpath('.'))

%% This is where the input would go...

% ---- CONDIUT ----
conIn = [];
% conIn.conduit_radius = 45.45; %39.313;
conIn.vh0 = 0;
% conIn.conduit_radius = 7.10879;
conIn.Q = 5e7;
conIn.conduit_radius = extrapVentRadius(conIn.Q);

% conIn.atmo = 'atm_ERAreanalysis_Shinmoedake2011_absWind.mat';
% conIn.atmo = 'atm_ERAreanalysis_Tungarahua2014_01_absWind.mat';
% conIn.atmo = 'atm_ERAreanalysis_Grimsvotn2011_01_absWind.mat';

% conIn.conduit_radius = 20.0034; %39.313;
% conIn.vh0 = 1700;
% conIn.Q = 9999750;

% conIn.atmo = 'atm_ERAreanalysis_Shinmoedake2011_absWind.mat';
% conIn.atmo = 'atm_ERAreanalysis_Tungarahua2014_01_absWind.mat';
conIn.atmo = 'atm_ERAreanalysis_Grimsvotn2011_01_absWind.mat';
% conIn.atmo = '~/Kahuna/data/tonga/atmoProfiles/HTHH_04_2022-01-15.mat';
% conIn.atmo = '~/Kahuna/data/gvolc-meghan/ERA5/Katla_avg_atmo_Oct12_1979-2005.mat';


% conIn.Q              = 1e8;
% conIn.Z0             = 6000;
% conIn.f0             = .001;
% conIn.Zw = 50;
% conIn.phi0           = 0.05;
% conIn.N0             = 1e9;
% conIn.dP = 0;
% conIn.conduit_radius = 37;
% conIn2.conduit_radius = 41.2;
% conIn2.Z0 = 5000;
% conIn2.Q = 1.0e8;
% conIn2.Q = 1e9;
% conIn2.Zw = 200;
% conIn2.conduit_radius = 44.6;
% conIn.Zw = 100;
% conIn.conduit_radius = 44.8;
conIn.proxy = true;
conIn.T = 1500;
conIn.n_0 = 0.02;
conIn.rho_melt = 2600;
conIn.vh0 = 850;
conIn2 = conIn;
% conIn2.n_ec = 0.1;
% conIn2.T = 900;
conIn2.Q = 5e7;
conIn2.conduit_radius = extrapVentRadius(conIn.Q);
conIn2.Zw = 70;

% ----  VAULT  ----
% blah blah

% ----  PLUME  ----
pS = [];
% pS.atmo = 'atm_ERAreanalysis_Grimsvotn2011_01_absWind.mat';
% pS.atmo = 'atm_ERAreanalysis_Shinmoedake2011_02_absWind.mat';
% rhob = 3.599205521505065; % Assumed bulk density
% pS.u_0 = 300;
% pS.r_0 = (1e7/(pi*pS.u_0*rhob)).^(1/2);
% pS.phi0 = 0;

% pI.useJetEntranceLength = false;
pS.D = 2.7;
pS.phiSz_max = 8;
pS.saScale   = 10;
pS.mwiModel  = 2;
pS.T_g_rng   = 50;
% pS.n_0 = .03;
% pS.phiSz_max = 10;
% pS.phiSz_min = 10; % Full pseudogas run
% pS.phi0      = 0;
%
pS2 = pS;
pS2.wind = false;
% pS2.phiFrag_cutoff = 3.43-1.46;
% pS2.useJetEntranceLength = false;
% pS2.alphaRTscale = 1;
% pS2.phiSz_min = -6;
% pS2.phiSz_max = 7;
% pS2.saScale   = 25;
% pS2.u_0 = 300;ed
% pS2.phi0 = 0;

%% DO THE THING 
tic
dat1 = runCoupledModel(conIn,pS);
toc;tic
dat2 = runCoupledModel(conIn2,pS); %,[],20);
toc

% --- OLD ROUTINE ---
% dat.conduitI  = getConduitSource;
% dat.conduitI  = getConduitSource(conIn);
% dat.conduitO  = Conduit_flow_with_nucleation(dat.conduitI);

% dat.plumeI  = con2plume(dat.conduitI,dat.conduitO);
% con2vault,vault2plume


% dat.plumeI.u0 = 300;
% dat.plumeI  = getPlumeSource(pS);
% dat.plumeO  = hmodel(dat.plumeI,true);
% dat2.plumeI  = getPlumeSource(pS2);
% dat2.plumeO  = hmodel(dat2.plumeI,true);

% ------------------

% Summary table
% singleRunSummary(dat)
% singleRunSummary(dat2)
%% Mass continuity check
% Conduit
% phi = dat.conduitO.porosity(end);
% % m_n_conduit = pi.*dat.conduitO.a(end).^2.*dat.conduitO.U(end).*(phi*dat.conduitO.rho_g(end) + (1-phi)*dat.conduitO.Par.rho_melt);
% m_n_conduit = pi.*dat.conduitO.a(end).^2.*dat.conduitO.U(end).*(dat.conduitO.rho_magma(end));
% % Plume
% m_0_plume = pi.*dat.plumeO.r(1).^2.*dat.plumeO.u(1).*dat.plumeO.rho_B(1);
% fprintf('Conduit exit mass flux:\t%.2f kg/s\n',m_n_conduit)
% fprintf('Plume source mass flux:\t%.2f kg/s\n',m_0_plume)

%% Plot
% plotCoupledPlume(dat)
% [axw,~] = plotMWIoutput_pap(dat1.wO,dat1.cI);
% plotMWIoutput(dat2.wO,dat2.cI,axw);
% plotMWIoutput([dat1 dat2],true,true)
% plotConduitPair([dat1 dat2])
plotMWIoutput([dat2],true,true)
plotPlumePair([dat1 dat2]);