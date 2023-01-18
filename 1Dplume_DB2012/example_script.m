clear all
close all


T0=1150; %exit temperature in K
n0=0.05; %exit gas content in wt.%
vh0=500; %vent altitude in m a.s.l.
u0=40; %exit velocity in m/s
r0=10;%vent radius in m


load('atmprofile.mat')
alt=atmprofile(:,2);%altitude in m
temp=atmprofile(:,3);%atmospheric temperature in K
hur=atmprofile(:,5)/100;%atmospheric relative humidity
wind=atmprofile(:,4);%atmospheric wind speed in m/s
press=atmprofile(:,1)*100;%patmospheric ressure in Pa

lambda=10^(-5);%condensation rate in s^-1
alpha=0.1%radial entrainment coefficient
beta=0.5%wind entrainment coefficient

u0_vec=zeros(50,1);
r0_vec=zeros(50,1);
for i=1:50
 %the last 3 inputs to the plume model are respectively: 1) the value of the
 %beta/alpha ratio used only for the 4th entrainment model; 2) the number of
 %entrainment model to be used (1 is for constant alpha and constant beta,
 %2 is variable alpha, 3 is variable alpha and variable beta, 4 is variable
 %alpha and constant beta/alpha ratio; see Aubry and Jellinek EPSL 2018;
 %and 3) the value of the normalization exponent as in Devenish et al (see
 %again our EPSL 2018 paper)
 
plumeSource = getPlumeSource('u0',u0,'r0',r0);
[pO] = hmodel(plumeSource);
% [htop rtop hb collapse mer] = hmodel(alt,temp,hur,wind,press,vh0,u0,T0,n0,r0,alpha,beta,lambda,NaN,1,1);
%htop = top height a.v.l.
%hb= NBL height a.v.l.
%collapse = boolean (0 for no collapse, 1 otherwise)
%mer = mass eruption rate in kg/s
%rtop=plume radius at the top in m
u0_vec(i)=u0;
r0_vec(i)=r0;

u0=u0*1.02;
r0=r0*1.09;
figure(1)
semilogx(pO.m_0,pO.hm,'ro')
hold on

semilogx(pO.m_0,pO.hb,'kv')


end
legend('Top height','NBL height')
xlabel('MER (kg/s)')
ylabel('Height (m a.v.l.)')

%
%
%
% for n0=0.02:0.01:0.07
% for beta=0.2:0.1:0.8
% [beta n0]
% for i=1:length(hlist)
%     load('atmprofile_eru.mat')
%     critmer=NaN(2,4,size(atmprofiles,3));
%
%         alt=squeeze(atmprofiles(:,2,82));
%         temp=squeeze(atmprofiles(:,3,82));
%         hur=squeeze(atmprofiles(:,5,82))/100;
%         wind=squeeze(atmprofiles(:,4,82));
%         press=squeeze(atmprofiles(:,1,82))*100;
%
%     [rho_B0 rho_aB0] = bdensity(alt,temp,hur,press,vh0,T0,n0);
%     ratio_ru=ricrit/(9.81*(rho_B0-rho_aB0)/rho_aB0);
% m0inilist=logspace(3,8,20);
% h0inilist=NaN(size(m0inilist));
% h0binilist=NaN(size(m0inilist));
% wfminilist=NaN(size(m0inilist));
% wfbinilist=NaN(size(m0inilist));
% %no condens
%
%
% for jj=1:length(m0inilist)
%
%     r0=((m0inilist(jj)/(pi*rho_B0))^2*ratio_ru)^(1/5);
%     u0=(m0inilist(jj)/(pi*rho_B0))/(r0^2);
%
%     [h0inilist(jj) rm wfminilist(jj) h0binilist(jj) wfbinilist(jj)] = hmodel(alt,temp,hur,wind,press,vh0,u0,T0,n0,r0,0.1,beta,lambda,NaN,1,1);
% end
%
%
%     h0inilist=h0inilist+vh0;
%     h0binilist=h0binilist*1.15+vh0;
%     mer_max_list(i)=interp1(h0inilist,m0inilist,hlist(i)*1000,'pchip');
%     mer_nbl_list(i)=interp1(h0binilist,m0inilist,hlist(i)*1000,'pchip');
%     wf_max_list(i)=interp1(h0inilist,wfminilist,hlist(i)*1000,'pchip');
%     wf_nbl_list(i)=interp1(h0binilist,wfbinilist,hlist(i)*1000,'pchip');
%
% end
%
%  figure(1)
%  subplot(2,2,1)
%  semilogy(hlist,mer_max_list,'ko','MarkerFaceColor','k')
%  hold on
%   subplot(2,2,2)
%  semilogy(hlist,mer_nbl_list,'ko','MarkerFaceColor','k')
%   hold on
%   subplot(2,2,3)
%  plot(hlist,wf_max_list,'ko','MarkerFaceColor','k')
%   hold on
%   subplot(2,2,4)
%  plot(hlist,wf_nbl_list,'ko','MarkerFaceColor','k')
%   hold on
%
%
% end
% end
%
% subplot(2,2,1)
% xlabel('Plume top height (km a.s.l.)')
% ylabel('Mass eruption rate at vent (kg/s)')
%
%
% subplot(2,2,2)
% xlabel('Plume spreading height (km a.s.l.)')
% ylabel('Mass eruption rate at vent (kg/s)')
%
% subplot(2,2,3)
% xlabel('Plume top height (km a.s.l.)')
% ylabel('water mass/ash mass at plume top')
%
%
% subplot(2,2,4)
% xlabel('Plume spreading height (km a.s.l.)')
% ylabel('water mass/ash mass at plume top')