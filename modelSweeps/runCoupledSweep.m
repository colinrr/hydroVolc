%% ============== Run coupledModel Sweep ===============
% Use conduitSweep output to run a set of couple model runs
%
% C Rowell, June 2021

clear all; close all;
% dDir     = '~/Kahuna/data/glaciovolc/';
% codeDir  = '~/code/research-projects/glaciovolc/glaciovolc-dev/';
dDir     = 'C:\Users\crowell\Kahuna\data\glaciovolc\';
codeDir  = 'C:\Users\crowell\Documents\GitHub\glaciovolc\glaciovolc-dev';
% dDir     = 'D:\Kahuna\data\glaciovolc\';
% codeDir  = 'C:\Users\crowell\Documents\GitHub\glaciovolc-dev';
addpath(genpath(codeDir))

conDir   = fullfile(dDir,'conduitSweeps/');
outDir   = fullfile(dDir,'coupledSweeps/');





outputName = 'coupledSweep_%s_MWIv%i_%s_n%i_%iq_%izw%s.mat'; % date, mwiVersion, atm_desc, numRuns, numQ, numZw

% Conduit run failure thresholds
ZfailScale  = 2;   % Z (depth) threshold in conduit radii
Mfailthresh = 0.95; % Mach number threshold
Pfailthresh = .05;  % Over-/underpressure threshold

% ------------ MODEL SWEEP PARAMS ------------
% Select subset of conduit runs (over surface water depth and MER)
% zwSubset = [0 20 50 100 200 400];
% q0Subset = 10.^[5.5:0.5:8.5];
% zwSubset = [0 100];
% q0Subset = 10.^[7 8];
% Full sweep
q0Subset = 10.^[5.5:0.1:8.9];
zwSubset = [];

% Run comparisons for these...eventually
% saScale  = [10 25 50];      % Particle roughness scale
% fragPar  = [0.05 0.1 0.2];  % Fragmentation energy partition

% ------------ COMMON MODEL PARAMS -----------
% ----- Conduit:

% High lat runs
% conIn.atmo = 'atm_ERAreanalysis_Grimsvotn2011_01_absWind.mat';
% conIn.vh0  = 1700;
qSweepFile = fullfile(conDir,'conduitV6_2021-06-16_Q_from_R-vs-Zw_Grimsvotn_V1.mat');
atm_desc   = 'hiLat';
plmIn.Tw0  = 274.15; % Surface water temperature           [K]
% desc       = ''; % Optional descriptor for output file name

% Low lat runs
% conIn.atmo = 'atm_ERAreanalysis_Tungarahua2014_01_absWind.mat';
% conIn.vh0  = 0;
% qSweepFile = fullfile(conDir,'conduitV6_2021-06-14_Q_from_R-vs-Zw_Tungarahua_V1.mat');
% atm_desc   = 'loLat';
% plmIn.Tw0  = 293.15; % Surface water temperature           [K]
% desc       = '';

% ----- MWI/Plume:
% Base run params
plmIn.mwiModel  = 2;      % Which magma-water-interaction model to use?
plmIn.D         = 2.9;    % PSD Power law exponent
plmIn.phiSz_max = 8;      % Max particle phi size (min particle dimension)
plmIn.saScale   = 10;     % Particle roughness scale
plmIn.T_g_rng   = 50;     % Range of temperatures for glass transition

% For RT entrainment runs (alphaRT_useLj)
plmIn.useJetEntranceLength = false;
plmIn.alphaRTscale = 0.5;
desc = '_alphaRT50_noLT';

% For noLd runs (noLd_freeD & noLd_cd) (getLd model 8 vs 9)
% plmIn.useDecompressLength = false;
% desc = '_noLd_cd';

% Control runs with no Ld and no Lj (noLd_noLj)
% plmIn.useDecompressLength = false;
% plmIn.useJetEntranceLength = false;
% desc = '_noLd_noLj';

% For noLd runs
% plmIn.useDecompressLength = false;
% desc = '_noLd';

% Control runs with no Ld and no Lj
% plmIn.useDecompressLength = false;
% plmIn.useJetEntranceLength = false;
% desc = '_noLd_noLj';

% Control runs with no Lj
% plmIn.useJetEntranceLength = false;
% desc = '_noLj';

% High Roughness runs
% plmIn.saScale = 25;
% desc = '_hiSA';

% fragPar runs
% plmIn.fragPar = 0.15;
% desc = '_hiFragPar15';

% hiD
% plmIn.D         = 3.2;
% desc = '_hiD';

%% ============== Do the thing ===============

% Load conduit sweep and get subset indices
nZw = length(zwSubset);
nQ  = length(q0Subset);

load(qSweepFile)

[zwChk,zwidx] = ismember(round(zwSubset),round(cQ.Zw));
[qChk,qidx]   = ismember(round(q0Subset),round([cI_base.Q]));
if ~isempty(zwSubset) && any(~zwChk)
    error('Zw subset elements are missing in conduit sweep')
end
if ~isempty(q0Subset) && any(~qChk)
    error('Q subset elements are missing in conduit sweep')
end
nRuns = nQ*length(cQ.Zw) + nZw*length([cI_base.Q]);

% Prep output?
Rmesh = repmat(cQ.R,[1 length(cQ.Zw)]);
qFmin = scatteredInterpolant(cQ.pf(cQ.success),Rmesh(cQ.success),cQ.Qmin(cQ.success),'linear','nearest');
qSetValid = zeros(nQ,length(cQ.Zw));
zSetValid = zeros(length(cQ.R),nZw);

% Run sweep
% Q subset along Zw first
ct = 0;
disp('Running sweeps along Zw...')
for ii=1:nQ
    qi = qidx(ii);
    for zi = 1:length(cQ.Zw)
        cI = cIq(qi,zi);
        pI = plmIn;
        
        % Check success, try interpolation
        if ~cQ.success(qi,zi)
            cI.Q = qFmin(cIq(qi,zi).pf,cIq(qi,zi).conduit_radius);
        end
        
        if ~isempty(cI.Q)
            ct = ct+1;
            tic
            qSet(ii,zi) = runCoupledModel(cI,pI);
            rTime = toc;
            
            [~,~,~,~,~,qSetValid(ii,zi)] = checkConduitResult(qSet(ii,zi).cO,cI.conduit_radius*ZfailScale,Mfailthresh,Pfailthresh);
        else
            rTime = 0;
            qSet(ii,zi).wO.failedPlume = true;
            qSet(ii,zi).pO.collapse = true;
            qSet(ii,zi).pO.hm = NaN;
            qSet(ii,zi).pO.hb  = NaN;
        end
        
        if and(ii==1,zi==1) % Allocate memory
            qSet(nQ,length(cQ.Zw)) = qSet(ii,zi);
        end
        
        % Spit out some output
        fprintf('  %i/%i:   runtime = %.0f s, Q = %.2e, Zw = %.0f, Cfail = %i, Wfail = %i, Clps = %i, Hm/Hm0 = %.2f\n',...
            ct,nRuns,rTime,cI.Q,cI.Zw,~qSetValid(ii,zi),qSet(ii,zi).wO.failedPlume,qSet(ii,zi).pO.collapse,qSet(ii,zi).pO.hm/qSet(ii,1).pO.hm)
    end
end
if isempty(q0Subset)
    qSet = [];
    qSetValid = [];
end

% Zw subset along Q next
disp('Running sweeps along Q...')
for jj = 1:nZw
    zi = zwidx(jj);
    for ri=1:length(cQ.R)
        cI = cIq(ri,zi);
        pI = plmIn;
        
        % Check success, try interpolation
        if ~cQ.success(ri,zi)
            cI.Q = qFmin(cIq(ri,zi).pf,cIq(ri,zi).conduit_radius);
        end
        if ~isempty(cI.Q)
            ct = ct+1;
            tic
            zSet(ri,jj) = runCoupledModel(cI,pI);
            rTime = toc;
            [~,~,~,~,~,zSetValid(ri,zi)] = checkConduitResult(zSet(ri,jj).cO,cI.conduit_radius*ZfailScale,Mfailthresh,Pfailthresh);
        else
            rTime = 0;
            zSet(ri,jj).wO.failedPlume = true;
            zSet(ri,jj).pO.collapse = true;
            zSet(ri,jj).pO.hm = NaN;
            zSet(ri,jj).pO.hb  = NaN;
        end
        
        if and(ri==1,zi==1) % Allocate memory
            zSet(length(cQ.R),nZw) = zSet(ri,jj);
        end
        
        % Spit out some output
        fprintf('  %i/%i:   runtime = %.0f s, Q = %.2e, Zw = %.0f, Cfail = %i, Wfail = %i, Clps = %i, Hm/Hm0 = %.2f\n',...
            ct,nRuns,rTime,cI.Q,cI.Zw,~zSetValid(ri,jj),zSet(ri,jj).wO.failedPlume,zSet(ri,jj).pO.collapse,zSet(ri,jj).pO.hm/zSet(ri,1).pO.hm)
    end
end
if isempty(zwSubset)
    zSet = [];
    zSetValid = [];
end

% Save output
oFile = sprintf(outputName,datestr(now,'yyyy-mm-dd'),plmIn.mwiModel,atm_desc,nRuns,nQ,nZw,desc);
oFile = fullfile(outDir,oFile);
fprintf('Saving:\n\t%s\n',oFile)
save(oFile,'qSet','zSet','zwSubset','q0Subset','qSweepFile','plmIn')

% Save summary output
[qA, oSumm] = getCoupledSweepArrays(oFile,true);