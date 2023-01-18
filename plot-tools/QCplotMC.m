% ####################################################
%               MONTE CARLO QC PLOTS
% ####################################################
clear all;  close all


dDir = '/Users/crrowell/Kahuna/data/gvolc-meghan/MonteCarloOutput/';

% MCfile = 'KatlaHydro_noLd_v1_2022-04-30_N2000.mat';
% MCfile = 'KatlaHydro_v1_2022-04-30_N2000.mat';
% MCfile = 'KatlaHydro_v2_2022-05-17_N2000.mat';
% MCfile = 'KatlaHydro_v3_2022-06-02_N1000.mat';
MCfile = 'KatlaHydro_v4_noLd_2022-06-23_N3500.mat';


% Plots - failure/success critera, summary data, randPars, raw data (indexed)?


addpath '/Users/crrowell/code/mybin/matlab/fig_scripts';
 %%
 
load(fullfile(dDir,MCfile)) %,'cI','pI','MC','randPars','summ')
load(fullfile(dDir,['outputSummary_' MCfile]))

%% Print summary info:
MCCf = fieldnames(MC.cI);
MCPf = fieldnames(MC.pI);

for ii=1:length(MCCf)
    ff = MCCf{ii};
    MCrep(ii).randPar = ff;
    MCrep(ii).dist    = MC.cI.(ff).dist;
    MCrep(ii).range_min    = min(MC.cI.(ff).range);
    MCrep(ii).range_max    = max(MC.cI.(ff).range);
end

disp('############# MONTE CARLO SIMULATION SUMMARY #############')
fprintf('File:\t%s\n',MCfile)
disp('Conduit Input:')
disp(cI)
disp('Plume Input:')
disp(pI)
disp('Random Params:')
disp(struct2table(MCrep))

Nsim = length(dat(:));
QClevel = zeros(length(dat),4);
QClevel(:,1) = all(qA.QCexist,2);
QClevel(:,2:4) = qA.QClevel;
QClevel = sum(QClevel,2);

fprintf('Sim. fail rate:      %i/%i (%.1f%%)\n',nansum(summ.modelFail),Nsim,nansum(summ.modelFail)/Nsim*100)
fprintf('No Output:           %i/%i (%.1f%%)\n',nansum(QClevel==0),Nsim,nansum(QClevel==0)/Nsim*100)
fprintf('Con - no vent sol.:  %i/%i (%.1f%%)\n',nansum(QClevel==1),Nsim,nansum(QClevel==1)/Nsim*100)
fprintf('MWI - failed breach: %i/%i (%.1f%%)\n',nansum(QClevel==2),Nsim,nansum(QClevel==2)./Nsim*100)
fprintf('Plume - error:       %i/%i (%.1f%%)\n\n',nansum(QClevel==3),Nsim,nansum(QClevel==3)./Nsim*100)

msg = cell(Nsim,1);
for ii=1:Nsim
    if ~isempty(summ.failMsg(ii).ME)
        msg{ii} = summ.failMsg(ii).ME.message;
    else
        msg{ii} = '';
    end
end
umsg = unique(msg);
nmsg = zeros(size(umsg));
disp('Detected errors:')
for jj=1:length(umsg)
    nmsg(jj) = sum(strcmp(msg,umsg(jj)));
    fprintf('  %4i : %s\n',nmsg(jj),umsg{jj})
end

%% QC plots

QCsymb = {'x','*','o','o','^'};

nr1 = 3;
nc1 = 2;
lw = 1.5;
dx = 0.1;
dy = 0.1;
ySz = [0.3 0.3 1];
ppads = [0.07 0.02 0.1 0.04];

% QCx1 = randPars.conduit_radius./extrapVentRadius(randPars.Q(:)); xl1 = 'a/a_0';
% QCx1 = log10(randPars.Q(:)); xl1 = 'log_{10}(Q)';
QCx1 = randPars.n_0(:) + randPars.n_ec(:); xl1 = 'n_w(z=0)';
% QCx1 = randPars.T(:);  xl1 = 'T (K)';
QCy1 = randPars.Zw(:)./randPars.conduit_radius; yl1 = 'Z_w/a';

QCx2 = qA.wO.z./randPars.Zw(:); xl2 = 'MWI: z/Z_e';
QCy2 = qA.wO.xv; yl2 = 'MWI: x_v(end)';


figure('Name','Model QC Checks')
tightSubplot(3,1,1,dx,dy,ppads,[],ySz)
imagesc(qA.QCexist')
colormap(flipud(jet))
set(gca,'YTick',1:6,'YTickLabel',{'cI','cO','pI','dO','wO','pO'})
title('Missing output')

tightSubplot(3,1,2,dx,dy,ppads,[],ySz)
imagesc(qA.QClevel')
set(gca,'YTick',1:3,'YTickLabel',{'Conduit','MWI','Plume'})
title('Simulation failures')

qcax = tightSubplot(nr1,nc1,5,dx,dy,ppads,[],ySz);
hold on
for ii=4:-1:0
    qcx = QCx1(QClevel==ii);
    qcy = QCy1(QClevel==ii);
    if ii>=3
        scatter(qcx,qcy,30,QCsymb{ii+1},'filled','LineWidth',lw)
    else
        scatter(qcx,qcy,30,QCsymb{ii+1},'LineWidth',lw)
    end
end
xlabel(xl1)
ylabel(yl1)
legend({'Success','Plume Fail','MWI Fail','Con. Fail','Empty'})

qcax = tightSubplot(nr1,nc1,6,dx,dy,ppads,[],ySz);
hold on
for ii=4:-1:0
    qcx = QCx2(QClevel==ii);
    qcy = QCy2(QClevel==ii);
    if ii>=3
        scatter(qcx,qcy,30,QCsymb{ii+1},'filled','LineWidth',lw)
    else
        scatter(qcx,qcy,30,QCsymb{ii+1},'LineWidth',lw)
    end
end
xlabel(xl2)
ylabel(yl2)
legend({'Success','Plume Fail','MWI Fail','Con. Fail','Empty'})