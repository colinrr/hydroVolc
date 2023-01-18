clear all;
dDir     = '~/Kahuna/data/glaciovolc/conduitSweeps';
% dDir     = 'D:\Kahuna\data\glaciovolc\conduitSweeps\';
fineFile = fullfile(dDir,'conduitV6_fineSweep_n64821_21-05-16_compressedV4_fineV2.mat');
% fineFile = fullfile(dDir,'conduitV6_fineSweep_n64821_21-05-16_compressedV4_fineV1.mat');

addpath('~/code/matlab/fig_scripts')
% addpath(genpath('C:\Users\crowell\Documents\Github\glaciovolc-dev\'))
% addpath('C:\Users\crowell\Documents\Github\pulseTracker-dev\pulseTracker\plot-tools')

ZfailScale = 2;   % Z (depth) threshold in conduit radii
Mfailthresh = 0.95; % Mach number threshold
Pfailthresh = .04;  % Overpressure threshold
Kthresh = 1.05;

% qSweepFile = fullfile(dDir,'conduitV6_2021-06-14_Q_from_R-vs-Zw_Tungarahua_V1.mat');
qSweepFile = fullfile(dDir,'conduitV6_2021-06-16_Q_from_R-vs-Zw_Grimsvotn_V1.mat');
%%
% load(fineFile)
load(qSweepFile)

%% Plot up the parameter space
fs = 12;
Q0 = cQ.Qmin(:,1);
logQ0 = log10(Q0);
dx = 0.05;
dy = 0.06;
ppads = [.04 .03 0.06 0.03];
nr = 3;
nc = 4;

%%
if 0==1
figure('position',[50 50 2000 1200])
% success, Rmax, Rmin, Ztop, M, K, rho
tightSubplot(nr,nc,1,dx,dy,ppads)
imagesc(Q0,cQ.Zw,double(cQ.success'))
shading flat
% view([0 0 1])
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'FontSize',fs)
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Success';

tightSubplot(nr,nc,2,dx,dy,ppads)
contourf(Q0,cQ.Zw,log10(cQ.Qmin'),logQ0)
% shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'log_{10}(Q_{min} (kg/s)) ';

tightSubplot(nr,nc,3,dx,dy,ppads)
contourf(Q0,cQ.Zw,log10(cQ.Qmax'),log10(cQ.Qmax(:,1)))
% shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'log_{10}(Q_{max} (kg/s)) ';

tightSubplot(nr,nc,4,dx,dy,ppads)
imagesc(Q0,cQ.Zw,cQ.Ztop')
shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Z_{min} (m)';
caxis([0 max(cQ.R(:))*2])

tightSubplot(nr,nc,5,dx,dy,ppads)
imagesc(Q0,cQ.Zw,cQ.U')
shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Velocity (m/s)';

tightSubplot(nr,nc,6,dx,dy,ppads)
imagesc(Q0,cQ.Zw,cQ.M')
shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Mach #';
hold on

Mcut = zeros(length(Q0)-1,1);
for ii=1:length(Mcut);Mcut(ii)= find(cQ.M(ii,:)>=Mfailthresh,1,'last');end
plot(Q0(1:end-1),smooth(cQ.Zw(Mcut),8),'LineWidth',2.5)

kax=tightSubplot(nr,nc,7,dx,dy,ppads);
K = cQ.pm./cQ.pf;
imagesc(Q0,cQ.Zw,log10(K)')
% view([0 0 1])
% axis tight
shading flat
% colormap jet
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'log_{10}(Overpressure)';
colormap(kax,redblue(100))
caxis(kax,(max(abs(log10(K(:))))*[-1 1]))
hold on

Kcut = zeros(length(Q0)-1,1);
for ii=1:length(Kcut);Kcut(ii)= find(K(ii,:)<=(Kthresh),1,'first');end
plot(Q0(1:end-1),smooth(cQ.Zw(Mcut),8),'LineWidth',2.5)
plot(Q0(1:end-1),smooth(cQ.Zw(Kcut),8),'LineWidth',2.5)

tightSubplot(nr,nc,8,dx,dy,ppads)
imagesc(Q0,cQ.Zw,cQ.rho_b')
shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Density (kg/m^3)';
hold on
plot(Q0(1:end-1),smooth(cQ.Zw(Mcut),8),'LineWidth',2.5)


tightSubplot(nr,nc,9,dx,dy,ppads)
imagesc(Q0,cQ.Zw,cQ.phi_g')
shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Vent gas vol. frac.';

% Fragmentation conditions
tightSubplot(nr,nc,10,dx,dy,ppads)
imagesc(Q0,cQ.Zw,cQ.Z_fr'/1e3)
shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Frag''n Depth (km)';

tightSubplot(nr,nc,11,dx,dy,ppads)
dPdt = cQ.dPdt_fr;
dPdt(~cQ.success) = NaN;
imagesc(Q0,cQ.Zw,dPdt'/1e6)
shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Frag''n dP/dt (MPa/s)';

tightSubplot(nr,nc,12,dx,dy,ppads)
imagesc(Q0,cQ.Zw,cQ.U_fr')
shading flat
colormap(gca,jet)
set(gca,'YDir','normal')
set(gca,'Xscale','log','XTick',10.^round(min(logQ0):max(logQ0)))
set(gca,'FontSize',fs)
xlabel('Q_0 (kg/s)')
ylabel('Z_w (m)')
cb = colorbar;
cb.Label.String = 'Frag''n Velocity (m/s)';

end

%% Show R successes as f'n of MER, PEF
if 0==1
figure
surf(MER,pf/1e6,cQ.Rmax','FaceAlpha',0)
cQgood = cQ.Rmax;
cQgood(~cQ.success) = NaN;
hold on
surf(MER,pf/1e6,cQgood','EdgeAlpha',0)
set(gca,'XScale','log')
xlabel('MER')
ylabel('P_a')
cb=colorbar;
colormap(jet(200))
end

%% Show evolution of radius search space for a fixed MER and several pressures

qi = 21; %1e7
% pi = [1 5:5:50];
pi = [1 2 6 10 20 50];

if 1==1
lalph = linspace(1,0.4,length(pi));

qmax = max([cQ.Qmax(qi,:)]);
rmin = min([cQ.Qmin(qi,:)]);
rrange = round(diff([rmin qmax]*0.1)*[-1 1] + [rmin qmax],1);
ri = rrange(1):0.1:rrange(2);

M = zeros(length(ri),length(pi));
K = M;
Z = M;
A = M;
Zf = M;
valid = false(size(M));

co = get(gca,'ColorOrder');
co = rgba2rgb(repmat(co(1,:),[length(lalph) 1]),lalph);
co = colormap(cool(length(lalph)));
if 0==1
    % cIt = cIf;
    % cIt.Q  = MER(qi);
    for ii=1:length(pi)
    %     cIt.pf = pf(pi);
        cIt = cI(qi,pi(ii));
         fprintf('P = %.3f MPa...\n',cIt.pf/1e6)

        % Run set of conduit runs
        for jj=1:length(ri)

            cIt.conduit_radius = ri(jj);
            cO = Conduit_flow_with_nucleation_V6(cIt);

            % Get valid checks
            Zfailthresh = ri(jj)*ZfailScale;
            [~,~,~,~,~,valid(jj,ii)] = checkConduitResult(cO,Zfailthresh,Mfailthresh,Pfailthresh);

            M(jj,ii) = cO.M(end);
            K(jj,ii) = cO.pm(end)./cO.Par.pf;
            Z(jj,ii) = cO.Z(end);
            Zf(jj,ii) = Zfailthresh;
            A(jj,ii) = cO.a(end)./cO.a(1);
        end

    end
    save(rSweepFile,'cI','cS','qi','pi','ri','M','K','Z','Zf','A','valid')
else
    load(rSweepFile)
end

end
%
%%

%% Plot Mach#, ZMin depth, Overpressure, max radius vs initial radius for up to 3 MER and 1 pf
if 1 ==1
ipick = 1:6; %5;
pi = pi(ipick);
M  = M(:,ipick);
K  = K(:,ipick);
Z  = Z(:,ipick);
Zf = Zf(:,ipick);
A  = A(:,ipick);
valid = valid(:,ipick);

lw1 = 1.5;
lw2 = 2.5;

figure('position',[50 50 1000 1100])
pa = pf(pi);

Zp = Z./Zf;
Zp(~valid) = NaN;
Z_rmax = cQ.Ztop(qi,pi)'./(cQ.Rmax(qi,pi).*2);

Kp = K;
Kp(~valid) = NaN;
K_rmax = (cQ.pm(qi,pi)./cQ.pf(qi,pi));

Mp = M; Mp(~valid) = NaN;
Ap = A; Ap(~valid) = NaN;
A_rmax = (cQ.amax(qi,pi)./cQ.Rmax(qi,pi));

win = cQ.success(qi,pi);

sp(1)=subplot(4,1,1);
plot(ri,ones(size(ri)),':k')
hold on
set(gca,'ColorOrder',co)
set(gca,'ColorOrderIndex',1)
plot(ri,Z./Zf,'--','LineWidth',lw1)
plot(ri,Zp,'.-','LineWidth',lw2)
scatter(cQ.Rmax(qi,pi(~win))',Z_rmax(~win)',60,co(~win,:),'LineWidth',lw1)
scatter(cQ.Rmax(qi,pi(win))',Z_rmax(win)',60,co(win,:),'filled')
ylabel('Min Depth (D/(Zf))')
ylim([0 10])
axis tight

sp(2) = subplot(4,1,2);
plot(ri,ones(size(ri))*(1+Pfailthresh),':k')
hold on
plot(ri,ones(size(ri))*(1-Pfailthresh),':k')
set(gca,'ColorOrder',co)
set(gca,'ColorOrderIndex',1)
plot(ri,K,'--','LineWidth',lw1)
set(gca,'ColorOrderIndex',1)
plot(ri,Kp,'.-','LineWidth',lw2)
scatter(cQ.Rmax(qi,pi(~win))',K_rmax(~win)',60,co(~win,:),'LineWidth',lw1)
scatter(cQ.Rmax(qi,pi(win))',K_rmax(win)',60,co(win,:),'filled')
ylabel('P_m/P_f')
axis tight

sp(3)=subplot(4,1,3);
plot(ri,ones(size(ri))*(Mfailthresh),':k')
hold on
set(gca,'ColorOrder',co)
plot(ri,M,'--','LineWidth',lw1)
set(gca,'ColorOrderIndex',1)
plot(ri,Mp,'.-','LineWidth',lw2)
scatter(cQ.Rmax(qi,pi(~win))',cQ.M(qi,pi(~win))',60,co(~win,:),'LineWidth',lw1)
scatter(cQ.Rmax(qi,pi(win))',cQ.M(qi,pi(win))',60,co(win,:),'filled')
ylabel('M')
axis tight

sp(4)=subplot(4,1,4);
plot(ri,ones(size(ri)),':k')
hold on
set(gca,'ColorOrder',co)
plot(ri,A,'--','LineWidth',lw1)
set(gca,'ColorOrderIndex',1)
plot(ri,Ap,'.-','LineWidth',lw2)
scatter(cQ.Rmax(qi,pi(~win))',A_rmax(~win)',60,co(~win,:),'LineWidth',lw1)
scatter(cQ.Rmax(qi,pi(win))',A_rmax(win)',60,co(win,:),'filled')
ylabel('a_v/a_0')
xlabel('a_0 (m)')
axis tight
linkaxes(sp,'x')

end
% set(gca,'ColorOrderIndex',1)
% plot([1;1]*cS.Rmin(iMER,pk)',[0;1]*[1 1],':')
% scatter(cS.Rmax(iMER,pk),cS.M(iMER,pk),50,'r')
% ylabel('M')
% legend({'7.5','8'})
% axis tight
% title(sprintf('MER = %.3e kg/s, P_f = %.3f MPa',MER(iMER),pf(pk)/1e6))
% 
% 
% plot(ric(:,iMER,pk),Ztc(:,iMER,pk)/1e3,'.-')
% hold on
% set(gca,'ColorOrderIndex',1)
% % plot([1;1]*cS.Rmax(iMER,pk)',[0;max(Ztc(:))]*[1 1],'--')
% set(gca,'ColorOrderIndex',1)
% % plot([1;1]*cS.Rmin(iMER,pk)',[0;max(Ztc(:))]*[1 1],':')
% scatter(cS.Rmax(iMER,pk),cS.Ztop(iMER,pk)/1e3,50,'r')
% ylabel('Z_{top} (km)')
% axis tight
% ylim([0 .1])
% 
% plot(ric(:,iMER,pk),pmc(:,iMER,pk)./pfc(:,iMER,pk),'.-');
% hold on
% set(gca,'ColorOrderIndex',1)
% plot([1;1]*cS.Rmax(iMER,pk)',[0;max(max(pmc(:,iMER,pk)./pfc(:,iMER,pk)))]*[1 1],'--')
% set(gca,'ColorOrderIndex',1)
% plot([1;1]*cS.Rmin(iMER,pk)',[0;max(max(pmc(:,iMER,pk)./pfc(:,iMER,pk)))]*[1 1],':')
% scatter(cS.Rmax(iMER,pk),cS.pm(iMER,pk)./pf(pk),50,'r')
% ylabel('K')
% xlabel('a (m)')
% axis tight
% linkaxes(sp,'x')

%% Run interpolations and check for solutions successes