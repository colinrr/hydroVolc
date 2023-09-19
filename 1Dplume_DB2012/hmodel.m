function varargout = hmodel(pI,verbose)
% pOut = hmodel(pI), [pOut,pIn] = hmodel(pS)]
% ---->
% pS = struct output of getPlumeSource
% verbose = true/false. true = output full z profiles. Otherwise just 5
%           basic output params
% pOut = output struct
% pIn  = optionally re-output pI (some calculated parameters will be
%          included)
%
% --> I/O changes, CR Mar 2021
% --> GSD implementation based on Girault ea 2014, Colucci ea 2014, CR Mar 2021

if nargin<2
    verbose = false;
end


alt         = pI.atmo(:,2);     % m a.s.l.
temp        = pI.atmo(:,3);     % K
hur         = pI.atmo(:,5)/100; % rel. hum
wind        = pI.atmo(:,4);     % m/s
press       = pI.atmo(:,1); %*100; % Pa
vent_height = pI.vh0;
u_0         = pI.u_0;
theta_0     = pI.T0;
n_0         = pI.n_0;
r_0         = pI.r_0;
alpha       = pI.alpha;
beta        = pI.beta;
omega       = pI.lambda;
baratio     = pI.baratio;
model       = pI.model;
nexp        = pI.nexp;

% Added for GSD - CR 2021
rho_m       = pI.rho_m; % Magma density. rho_s is now bulk pyroclast density
ri          = pI.Rgsd;  % Particle radii (m)
rhoi        = pI.rhoi;  % Bulk density of pyroclasts at each grain size
ni          = pI.ni;    % Fraction of pyroclast mass in bubble gas at each grain size
nsi_0        = pI.nsi;    % Fraction of total pyroclast mass at each grain size
numGS       = length(ri); % Number of grain size bins;

% Exported to I/O - CR 2021
% rho_s       = rho_m; % Null porosity case
rho_s       = (sum(pI.nsi./pI.rhoi))^-1; % Bulk density over GSD
C_s         = pI.C_s;

% Water props and MWI params - CR 2021
xv_0        = pI.xv_0; % Water dryness fraction


alt     = alt(~isnan(alt));
temp    = temp(~isnan(alt));
hur     = hur(~isnan(alt));
press   = press(~isnan(alt));
if pI.wind
    wind    = wind(~isnan(alt));
else
    wind = zeros(size(~wind(~isnan(alt))));
end
% whos

%constants
g          = 9.81;  % gravitational acceleration (m/s^2)
C_d        = 998;   % specific heat of air at constant pressure (J/kg/K)
C_v        = 1952;  %2000 % specific heat of volcanic gas (water) at constant pressure (J/kg/K) (@~500 K)
C_l        = 4190;  % specific heat of liquid water (J/kg/K)
% C_s        = 1250;  %850 % specific heat of solid pyroclasts (J/kg/K)
R_d        = 287;   % gas constant of dry air (J/kg/K) 
R_v        = 461;   % gas constant of volcanic gas (water) (J/kg/K) 
L          = 2.257e6;  % latent heat of vaporization (J kg^-1)
rho_l      = 1000;  % density liquid water (kg/m^3)
% rho_s      = 2350;  %1000 % density of solid pyroclasts (kg/m^3)
mu_a       = 1.8e-5; % Air dynamic viscosity (small variation w/ T,P) (Pa s)
xv_0       = pI.xv_0; % Steam/water dryness fraction

z0         = 0; % start calculation at vent height


% omega      = 0; % 0, no condensation, 1/(2.8*3600) moderate, 1/(1.7*60) rapid
% % omega matters even in a dry environment because water is a priori present
% % in plume
% r_h = 0; % relative humidity

% Vmax       = 10;

% vent_height       = 1500; % vent height (m)
Meteo_Height      = double(alt)- vent_height;  % Set Z=0 at vent height
Meteo_Temperature = double(temp);
Meteo_Wind        = double(wind);
Meteo_Pressure    = double(press);
Meteo_Humidity    = double(hur);
P_0               = 10.^interp1(Meteo_Height,log10(Meteo_Pressure),0,'pchip','extrap'); % ambient pressure at vent height
theta_a_0 = interp1(Meteo_Height,Meteo_Temperature,0,'pchip','extrap');
hur_a_0= interp1(Meteo_Height,Meteo_Humidity,0,'pchip','extrap');

% Water props
% bulk water density
if xv_0 == 1
    rho_v0 = density(P_0,theta_0);
else
    rho_v0 = densityTX(theta_0,xv_0);
end
% rho_v0  = P_0/R_v/theta_0; %P_0/(R_v*theta_0)

phi_s0  = (1 - n_0)*rho_v0/(n_0*rho_s+(1-n_0)*rho_v0);

rho_B0  = (n_0/rho_v0 + (1-n_0)/rho_s)^(-1);
eps_0   = R_d/R_v;
es_0    = 100.*6.112*exp(17.67*(theta_a_0-273.15)./(theta_a_0 +243.5-273.15)); % saturation vapor pressure for water
w_s_0   = 1/eps_0 * es_0/(P_0-es_0); % w_a = r_h*w_s, mass mixing ratio of water vapor to dry air at saturation
w_a_0   = hur_a_0*w_s_0; % mass mixing ratio of water vapor to dry air
rho_aB0 = P_0/(R_v*theta_a_0)*(1+w_a_0)/(w_a_0+eps_0);






phi_v0   = n_0*rho_s/(n_0*rho_s+(1-n_0)*rho_v0);

m_0      = rho_B0*u_0*r_0^2;
v_0      = pi*u_0*r_0^2;
b_0      = pi*g*((rho_B0-rho_aB0)/rho_aB0)*u_0*r_0^2;

x_0      = 0;
z_0      = z0;
m_d0     = 0;
m_v0     = n_0*(xv_0)*m_0;
m_l0     = n_0*(1-xv_0)*m_0;
m_s0     = rho_s*phi_s0*m_0/rho_B0;
m_si0    = m_s0.*nsi_0;                % Particle mass at each grain size
% m_0      = m_d0 + m_v0 + m_l0 + m_s0;
psi_0    = m_0*u_0;
angle_0  = 0.5*pi;
C_sb0    = C_v.*sum(ni.*nsi_0) + C_s*sum((1-ni).*nsi_0); % Initial bulk pyroclast heat capacity
C_B0     = (m_d0*C_d + m_v0*C_v + m_l0*C_l + m_s0*C_sb0)/m_0;
Q_0      = m_0*C_B0*theta_0;

Ri0=g*((rho_aB0-rho_B0)/rho_aB0)*r_0/(u_0^2);
De=2*r_0;
Lme=pi^0.25*r_0*abs(Ri0)^(-0.5);
wse=10.^interp1(Meteo_Height,log10(Meteo_Wind),0,'pchip','extrap')/u_0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   solve system of ODE's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('RelTol',1e-7,'AbsTol',1e-6,'Events',@(z,y) stopHybrid(y,numGS,u_0));
%     %error tolerance, and stop at velocity=0m/s (see stopPlume.m file)
Sspan = [z_0 100000]; % solve from z=0m to z=40000m
IC = [x_0 z_0 m_d0 m_v0 m_l0 m_si0' psi_0 angle_0 Q_0]; % initial conditions, add P_0 if hydrostatic approx in ode

odeSolution = ode15s(@(s,y) odeHybrid(y, ...
                         g, ...
                         C_d, ...
                         C_v, ...
                         C_l, ...
                         C_s, ...
                         R_d, ...
                         R_v, ...
                         L,   ...
                         model,      ...
                         nexp,      ...
                         alpha,      ...
                         beta,       ...
                         baratio,       ...
                         omega,      ...
                         mu_a,      ...
                         ri,        ...
                         rhoi,      ...  % m_s0,       ...
                         ni,        ...
                         De,       ...
                         Lme,       ...
                         u_0,       ...
                         vent_height,       ...
                         Meteo_Humidity,        ...
                         Meteo_Height, ...
                         Meteo_Temperature, ...
                         Meteo_Wind, ...
                         Meteo_Pressure),Sspan,IC,options);
                     
                     
S = odeSolution.x';
Y = odeSolution.y';

 if ~isreal(sum(sum(Y)))
        hm=NaN;
 rm=NaN;
 else
x         = Y(:,1);     % x-coordinate, plume centerline
z         = Y(:,2);     % z-coord
m_d       = Y(:,3);     % dry air mass flux?
m_v       = Y(:,4);     % H20 vapour mass flux
m_l       = Y(:,5);     % H20 liquid mass flux
m_si        = Y(:,6:6+numGS-1); % Solids mass flux at each grain size
m_s       = sum(m_si,2);       % Solids mass flux
psi       = Y(:,6+numGS);     % Momentum flux
angle     = Y(:,7+numGS);     % Centerline angle
Q         = Y(:,8+numGS);     % Buoyancy flux
%P         = Y(:,9);

nsi        = m_si./m_s;

eps       = R_d/R_v;
m         = m_d + m_v + m_l + m_s;  % Total mass flux
u         = psi./m;
% display(m)
C_sb      = C_v.*sum(ni'.*nsi,2) + C_s*sum((1-ni').*nsi,2);
C_B       = (m_d.*C_d + m_v.*C_v + m_l.*C_l + m_s.*C_sb)./m;
% display(C_B)
theta     = (1./C_B).*(Q./m);
rho_B     = (10.^interp1(Meteo_Height,log10(Meteo_Pressure),z,'pchip','extrap')./R_v./theta).*m./(m_v+eps.*m_d);
r         = (m./(rho_B.*u)).^(1/2);



    [hm ind]=max(z);
     rm=r(ind);




%buoyancy
eps       = R_d/R_v;
es = 100.*6.112*exp(17.67*(interp1(Meteo_Height,Meteo_Temperature,z,'spline','extrap')-273.15)./(interp1(Meteo_Height,Meteo_Temperature,z,'spline','extrap') +243.5-273.15)); % saturation vapor pressure for water
w_s = 1/eps * es./(10.^interp1(Meteo_Height,log10(Meteo_Pressure),z,'spline','extrap')-es); % w_a = r_h*w_s, mass mixing ratio of water vapor to dry air at satuartion
w_a = interp1(Meteo_Height,Meteo_Humidity,z,'spline','extrap').*w_s; % mass mixing ratio of water vapor to dry air



collapse=0;

rho_aB    = 10.^interp1(Meteo_Height,log10(Meteo_Pressure),z,'spline','extrap')./(R_v*interp1(Meteo_Height,Meteo_Temperature,z,'spline','extrap')).*(1+w_a)./(w_a+eps);
drho=rho_B-rho_aB;
sdrho=sign(drho);
rhoinv=[0;diff(sdrho)];
indinv=find(rhoinv~=0);
if numel(indinv)==0
    collapse=1;
    hb=NaN;
   
elseif drho(1)<0 & numel(indinv)==1
     hb=interp1([drho(indinv(1)-1) drho(indinv(1))],[z(indinv(1)-1) z(indinv(1))],0);
    
     
else
    hb=interp1([drho(indinv(2)-1) drho(indinv(2))],[z(indinv(2)-1) z(indinv(2))],0);
  
end

m_0      = pi*m_0 ; % CR 2021....wtf?
   
% ---- new output setup CR 2021 -----
Pout.hm = hm;               % plume top height  [m]
Pout.rm = rm;               % plume radius at top [m]
Pout.hb = hb;               % neutral bouyancy height [m]
Pout.collapse = collapse;   % column collapse flag
Pout.m_0 = m_0;             % source mass eruption rate [kg/s]

if isempty(pI.rho_B0)
    pI.rho_B0 = rho_B0;
end
if isempty(pI.rho_g0)
    pI.rho_g0 = rho_v0;
end

if verbose
    ztropo = findTPheight(Meteo_Height/1e3,Meteo_Temperature);

     if max(real(z))<=max(Meteo_Height)
        Pout.atmo.theta_a = interp1(Meteo_Height,Meteo_Temperature,real(z),'pchip','extrap');
    else
        Pout.atmo.theta_a = interp1([Meteo_Height;50000-vent_height;85000-vent_height],[Meteo_Temperature;273;173],real(z),'pchip','extrap');
    end

    
    Pout.atmo.es    = es;
    Pout.atmo.w_s   = w_s;
    Pout.atmo.w_a   = w_a;
    Pout.atmo.rho   = rho_aB;
    Pout.atmo.ztropo = ztropo*1e3; %findTPheight(z/1e3,Pout.atmo.theta_a); %z(find(gradient(Pout.atmo.theta_a,z/1e3)>-2,1,'first'));
    Pout.atmo.zDatum = vent_height;
    Pout.P_0    = P_0;
    Pout.angle  = angle;
    Pout.C_B    = C_B;
    Pout.drho   = drho;
%     Pout.es     = es;
    Pout.m      = m*pi;
    Pout.m_d    = m_d*pi;
    Pout.m_l    = m_l*pi;
    Pout.m_v    = m_v*pi;
    Pout.m_si   = m_si*pi;
    Pout.nsi     = nsi;
    Pout.m_s    = m_s*pi;
    Pout.psi    = psi;
    Pout.Q      = Q;
    Pout.r      = r;
%     Pout.rho_aB = rho_aB;   
    Pout.rho_B  = rho_B;
    Pout.rhoinv = rhoinv;
    Pout.S      = S;
    Pout.sdrho  = sdrho;
    Pout.theta  = theta;
    Pout.u      = u;
%     Pout.w_a    = w_a;
%     Pout.w_s    = w_s;
    Pout.z      = z;
%     Pout.drho
end

varargout{1} = Pout;
if nargout==2
    varargout{2} = pI;
end
% -----------------------------------
    end
    
   
   
    



% figure(1)
% subplot(2,2,4), plot(rho_aB,z.*1e-3,'b')
