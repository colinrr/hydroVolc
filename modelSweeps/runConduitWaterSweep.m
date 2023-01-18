%% ============== Run coupledModel Sweep ===============
% Use conduitSweep output to run a set of couple model runs
%
% C Rowell, June 2021

clear all; close all;

% dDir     = '~/Kahuna/data/glaciovolc/';
% codeDir  = '~/code/research-projects/glaciovolc/glaciovolc-dev/';
% dDir     = 'C:\Users\crowell\Kahuna\data\glaciovolc\';
% codeDir  = 'C:\Users\crowell\Documents\GitHub\glaciovolc\glaciovolc-dev';
dDir     = 'D:\Kahuna\data\glaciovolc\';
codeDir  = 'C:\Users\crowell\Documents\GitHub\glaciovolc-dev';
addpath(genpath(codeDir))

outDir   = fullfile(dDir,'conduitWaterSweeps/');
outputName = 'conWaterSweep_%s_MWIv%i_PROXv%i_%s_n%i_%iq_%izw%s.mat'; % date, mwiVersion, ventProxyVersion, atm_desc, numRuns, numQ, numZw
desc = '_testRun';

qi = 19; % Which MER to use? 19 = 1e7, 29 = 1e8
cS.Q0 = 1.9952e+07;
cS.conduit_radius = 25.5849;

cS.vh0  = 1700;
cS.atmo = 'atm_ERAreanalysis_Grimsvotn2011_01_absWind.mat';
% cS.atmo = 'atm_ERAreanalysis_Tungarahua2014_01_absWind.mat';
atm_desc = 'hiLat';
beta_crit = 0.95;

Ze   = [0:5:300]';
n_ec = [0 0.02 0.05 0.1 0.15]; 

% Base run params - Reference Scenario
plmIn.mwiModel  = 2;      % Which magma-water-interaction model to use?
plmIn.D         = 2.9;    % PSD Power law exponent
plmIn.phiSz_max = 8;      % Max particle phi size (min particle dimension)
plmIn.saScale   = 10;     % Particle roughness scale
plmIn.T_g_rng   = 50;     % Range of temperatures for glass transition
plmIn.Tw0       = 274.15; % Surface water temperature           [K]

%%
nRuns = length(Ze)*length(n_ec);
nZw = length(Ze);

ct = 0;
disp('Running sweep...')
for zi=length(Ze):-1:1
    
    for ni = 1:length(n_ec) %:-1:1
        ct = ct+1;
        cI = cS;
        pI = plmIn;
        
        cI.Zw = Ze(zi);
        cI.n_ec = n_ec(ni);

        % Spit out some output
        if isfield(dat(zi,ni).wO,'n_l')
            w_liq = dat(zi,ni).wO.n_l(end)/(dat(zi,ni).wO.n_l(end) + dat(zi,ni).wO.n_v(end));
            Zf    = dat(zi,ni).wO.z(end)/dat(zi,ni).cI.Zw;
        else
            w_liq = 0;
            Zf = NaN;
        end
        
        % Spit out some output
        w_liq = dat(zi,ni).wO.n_l(end)/(dat(zi,ni).wO.n_l(end) + dat(zi,ni).wO.n_v(end));
        Zf    = dat(zi,ni).wO.z(end)/dat(zi,ni).cI.Zw;
        fprintf('  %i/%i:   runtime = %.0f s, Q = %.2e, Zw = %.0f, n_ec = %.2f, Choke = %i, Wfail = %i, %%w_liq = %.2f, %%Zw_f = %.2f, Clps = %i, Hm/Hm0 = %.2f\n',...
            ct,nRuns,rTime,cI.Q0,cI.Zw,cI.n_ec,dat(zi,ni).cO.Par.choke,dat(zi,ni).wO.failedPlume,w_liq,Zf,dat(zi,ni).pO.collapse,dat(zi,ni).pO.hm/dat(zi,1).pO.hm)
        
    end
end

% Save output
oFile = sprintf(outputName,datestr(now,'yyyy-mm-dd'),plmIn.mwiModel,2,atm_desc,nRuns,1,nZw,desc);
oFile = fullfile(outDir,oFile);
fprintf('Saving:\n\t%s\n',oFile)
save(oFile,'dat','cS','plmIn','beta_crit','Ze','n_ec')
 
% Save summary output
[qA, oSumm] = getProxySweepArrays(oFile,true);

%% Fig params

%% Make fig