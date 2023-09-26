

function [output] = Conduit_flow_with_nucleation_V7(Input)
% By Sahand Hajimirza
% Last update: Apr 22, 2021

% This code calculates bubble nucleation and growth during magma ascent in
% a cylindrical conduit with a constant cross sectional area. 
% The model description and equations are given in 
% Hajimirza et al "The shape of volcanic conduits inferred from bubble size
% distributions" EarthArxiv 2020

% Updates to V6 introduced by Colin Rowell
% --> Included optional composition input to tweak viscosity model
% --> Included workaround/fix for cases with initial non-zero bubble number density, Sep 2023
% --> Minor (mostly I/O) modifications to incorporate getConduitSource function, Apr/May 2021
% --> V2/V6 versioning uncertain here?

% Updates to V2:
% Flaring has been added
% Flaring is calculated using Mastin 2002 formulation with speed of sound

% Updates to V1:
% Speed of sound in the gas phase has been corrected to
% account for isothermal flow.

% Sound velocity after fragmentation has been changed from in pure gas to in
% the mixture (psudo-gas). 

global fw_interpolant psat_interpolant 

% Check for proper input structure
assert(all(isfield(Input, {'atmo','composition','conduit_radius','dP','f0',...
    'Mfailthresh','N0','pf','Pfailthresh','phi0','phi_frag','proxy','Q',...
    'rho_melt','ST_coeff','T','theta','vh0','Z0','ZfailScale','Zw'})),...
    'Conduit model input is missing fields - check getConduitSource')

%==========================================================================
% Physical properties and constants
%==========================================================================
Par.Mw          = 18.0528e-3;                   % Molar Mass of water [Kg/mol]
Par.KB          = 1.38e-23;                     % Boltzman constant [J/K]
Par.AV          = 6.022e23;                     % Avogadro number [1/mol]
Par.Rgas        = 8.314;                        % Gas constant [J/(K.mol)]
% Par.rho_melt    = 2400;                         % Melt density [Kg/m^3]
Par.g           = 9.81;                         % Gravity
rho_rock        = 2400;                         % Rock density
% Par.f0          = .0025;                        % Friction coefficient (Mastin 2000)
Par.K_melt      = 224e8;                        % Bulk Modulus of melt

% --- CR 2021 input edits ---
Par.rho_melt = Input.rho_melt;
Par.pf       = Input.pf;
Par.f0       = Input.f0;
frag         = 0;
choke        = 0;
Zf           = NaN;
Zc           = NaN;
% ---------------------------

%==========================================================================
% Initial Conditions
%==========================================================================
Par.Z0          = Input.Z0;                                 % Initial depth
% Par.Pinitial    = rho_rock * Par.g * Par.Z0 + Input.dP;     % Initial pressure = Lithostatic + overpressure
% Par.pf          = 1e5;                                      % Final pressure = atmpspheric pressure
Par.Pinitial    = Par.pf + rho_rock * Par.g * Par.Z0 + Input.dP; % Initial pressure = Lithostatic + overpressure
Par.T           = Input.T; % + 273.15;                         % Temperature (CR: switched to K)
Par.Xc          = 0;                                        % Mole fraction of CO2
Par.Q           = Input.Q;                                  % Mass discharge rate
Par.C0          = solubility(Par.Pinitial,Par.Xc,Par.T);    % Initial water concentration
Par.mu          = viscosity(Input.composition,Par.T,Par.C0);
Par.ST_coeff    = Input.ST_coeff;                           % Heterogeneous nucleation factor
Par.phi_frag    = Input.phi_frag;
Par.P_frag      = 1e5;
Par.a           = Input.conduit_radius;                     % Conduit radius
Par.composition = Input.composition;
         
Par.Mc = .97;   %  Critical Mach number

% Interpolation for fugacity coefficient as a function of pressure at 850C.
% It can be written as a function of pressure and temperature if calculation at other
% temperatures is needed.
load('fw_interpolant.mat','fw_interpolant_850C')
fw_interpolant = fw_interpolant_850C;

% Interpolation for finding saturation pressure as a function of water
% concentration at 850 C. A general function as function of temperature can
% be written.
load('find_psat_interpolant.mat')
psat_interpolant = psat_interpolant;




%===========================================================================
Par.ND  = 1;            % Initial Threshold after which bubble growth is calculated
M0      = Input.N0;     % Number density of initial exsolved volatiles

% If there is no initial exsolved volatiles, the model calculates
% nucleation rate until number of bubbles reaches the minimum threshold,
% After which the code calculates bubble nucleation and growth together.

% If there is initial exsolved volatiles, the code directly goes into
% bubble nucleation and growth calculations. 

if M0 < Par.ND

    option1=odeset('AbsTol',1e-5,'RelTol',1e-5,...
        'Events',@EventsFcn1, 'MaxStep',10);
    y0 = [M0 Par.Pinitial Par.Z0]';     % Initial BND
    [t1,y1,~,~,~] = ode15s(@before_nucleation,[0 1e10],y0,option1,Par);

    [~,par1]=before_nucleation(t1,y1',Par);
    te1 = par1.t(end);
    
    
    % Preparing initial conditions for bubble nucleation and growth
    % function

    M0 = par1.M0(end);   
    r = 2 * par1.rc(end);
    M1 = M0 * r;
    M2 = M0 * r^2;
    M3 = M0 * r^3;
    pg = par1.pm(end) + 2*par1.sigma(end)/par1.rc(end);
    mg = (pg*4/3*pi*r^3*Par.Mw) / (Par.Rgas*Par.T)*M0;
    Cm = Par.C0 - M0/Par.rho_melt * mg;
    pm = par1.pm(end);
    z = par1.Z(end);

else
    
    % Preparing initial conditions for bubble nucleation and growth
    % function
    M0 = Input.N0;   
    M3 = 3/4/pi * Input.phi0/(1-Input.phi0);
    r = (M3/M0)^(1/3);
    M1 = M0 * r;
    M2 = M0 * r^2;

    pressure.pm = Par.Pinitial;
    pressure.psat = pressure.pm;
    pressure.pb = findPb(pressure.pm,pressure.psat,Par.T,fw_interpolant);
    ST = SurfaceTension(pressure,Par.T);

    pg = Par.Pinitial + 2*ST/r;
    rhog = EoS_H2O_2(pg,Par.T);
    mg = rhog * 4/3*pi*M3;

    Par.C0 = solubility(Par.Pinitial+Input.dP,Par.Xc,Par.T) + mg/Par.rho_melt;

    pm = Par.Pinitial;
    z = Par.Z0;

    te1 = 0;
    
    
end


%==========================================================================
% Before fragmentation (bubble nucleation and growth)
%==========================================================================

y0 = [mg,M0,M1,M2,M3,pm,z]';     % Initial conditions

options=odeset('AbsTol',[1e-38,1e-3,1e-12,1e-20,1e-28,1e-4,1e-5],...
               'RelTol',1e-5,...
               'Events',@EventsFcn2);
[t2,y2] = ode15s(@before_fragmentation,[te1 1e10],y0,options,Par);
[~,par2] = before_fragmentation(t2,y2',Par);

%==========================================================================
% After fragmentation (Gas with dispersed pyroclasts)
%==========================================================================
if par2.Z(end) > 100 && par2.pm(end)/1e6 > .5
   
    y0 = [par2.mg(end),par2.M0(end),par2.M1(end),par2.M2(end),par2.M3(end),...
        par2.pm(end),par2.Z(end),par2.pg(end), par2.a(end)];
    options=odeset('AbsTol',[1e-25,1e-3,1e-12,1e-20,1e-28,1e-2,1e-5,1e-2, 1e-2],...
               'RelTol',1e-5,...
               'Events',@EventsFcn3);
    [t3,y3] = ode15s(@after_fragmentation,[par2.t(end) 1e10],y0,options,Par);
    [~,par3] = after_fragmentation(t3,y3',Par);
    
    if par3.M(end) >= Par.Mc && par3.Z(end) > 10
        y0 = [par3.mg(end),par3.M0(end),par3.M1(end),par3.M2(end),par3.M3(end),...
        par3.pm(end),par3.Z(end),par3.pg(end), par3.a(end)];
        options=odeset('AbsTol',[1e-25,1e-3,1e-12,1e-20,1e-28,1e-2,1e-6,1e-2, 1e-2],...
               'RelTol',1e-5,...
               'Events',@EventsFcn4);
        [t4,y4] = ode15s(@flaring,[par3.t(end) 1e10],y0,options,Par);
        [~,par4] = flaring(t4,y4',Par);
        
        % preparing output
        parnames = fieldnames(par3);
        for i = 1:length(parnames)
       
            par3.(char(parnames(i))) =...
            [par3.(char(parnames(i))); par4.(char(parnames(i)))];
        end
        
        % --- CR edit: record choking flag/depth 
        choke = 1;
        Zc    = par3.Z(end);
    end

%==========================================================================

% Preparing output 
    parnames = fieldnames(par2);

    for i = 1:length(parnames)
       
        if exist('par1','var') % Workaround for case with initial BND, in which par1 (pre-nucleation) is skipped
            output.(char(parnames(i))) =...
                [par1.(char(parnames(i))); par2.(char(parnames(i)));...
                par3.(char(parnames(i)))];
        else
            output.(char(parnames(i))) =...
                [par2.(char(parnames(i)));...
                par3.(char(parnames(i)))];
        end
    end

    % CR edit: record fragmentation flag/depth
    frag = 1; 
    Zf   = par2.Z(end);
    %----
else
    
    parnames = fieldnames(par2);

    for i = 1:length(parnames)

        if exist('par1','var') % Workaround for case with initial BND, in which par1 (pre-nucleation) is skipped
            output.(char(parnames(i))) =...
            [par1.(char(parnames(i))); par2.(char(parnames(i)))];
        else
            output.(char(parnames(i))) =...
                par2.(char(parnames(i)));
            
        end
 
    end

    % CR edit: fragmentation flag/depth
%     frag = 0;
%     Zf   = NaN;
    %----
end    
output.Csol = solubility(output.pm,Par.Xc,Par.T);

% Added by CR (2021) for a few additional output params
output.Par    = Par;
output.Par.Zf = Zf;     % Fragmentation depth
output.Par.Zc = Zc;     % Choking depth
output.Par.frag = frag; % Fragmentation flag
output.Par.choke = choke; % Choking flag
% -----
end



function [dydt,out]=before_nucleation(t,y,Par)

global fw_interpolant 

    out.t = t;
    out.M0 = [y(1,:)]';
    out.pm = [y(2,:)]';
    out.Z = [y(3,:)]';
    


    %=============================================
    % Calcuating bubble nucleation rate
    out.Cm = Par.C0 * ones(size(t));
    out.D = DH2Orhyolite(out.Cm,out.pm,Par.T);    %Diffusivity of volatile in melt [m^2/s]
    out.psat = Par.Pinitial * ones(size(t));
    out.pb = findPb(out.pm,out.psat,Par.T,fw_interpolant); 
    out.sigma = SurfaceTension(out,Par.T); 
    
    out.rc = 2 * out.sigma ./ (out.pb - out.pm);    % Critical bubble radius
    n0 = out.Cm * Par.rho_melt * Par.AV / Par.Mw;
    a0 = n0.^(-1/3);
    Vw = Molecular_VH2O(out.pm,Par.T);
    J0 = (2*Vw .* n0.^2 .* out.D) ./ (a0) .* sqrt(out.sigma./(Par.KB*Par.T));
    
    out.Wcl = (16*pi * out.sigma.^3) ./ (3 * (out.pb-out.pm).^2); % Nucleation work
    
    out.J = J0 .* exp(-out.Wcl./(Par.KB*Par.T) * Par.ST_coeff);    % Nucleation rate
    %=============================================
    
    
    out.mu = viscosity(Par.composition,Par.T,Par.C0) * ones(size(t));       % Viscosity
    out.rho_magma = Par.rho_melt * ones(size(t));           % Density

    
    out.a = Par.a * ones(size(t));      % Conduit radius
    out.dadz = zeros(size(t));          % Conduit radius derivative in respect to z
    
    
    K_magma = Par.K_melt;                       % Bulk modulus of magma
    out.C = (K_magma ./ out.rho_magma).^(1/2);  % Sound velocity
    U_ave = Par.Q./(out.rho_magma*pi.*out.a.^2);
    out.U = U_ave;
    out.M = out.U ./ out.C;     % Mach number

    
    %=============================================
    % Calcuating decompression
    % Friction factor, Mastin 2000 eq. 18
    fric = 16*out.mu ./ (out.rho_magma.*U_ave*2.*out.a) + Par.f0;  

    % Decompression, Mastin 2000 eq. 7
    dPdz = (out.rho_magma * Par.g ...
     + out.rho_magma .* U_ave.^2 .* fric./out.a ...
     - 2 * out.rho_magma .* U_ave.^2 ./ out.a .* out.dadz)./(1-out.M.^2);

    
    out.dPdt = dPdz .* U_ave;
    
    dydt = [out.J; -out.dPdt; -out.U];
    
    
    
    %=============================================
    % Output variables, must be consistent with the functions with bubble
    % growth
    out.cin = NaN * ones(size(t));
    out.pg = NaN * ones(size(t));
    out.dC_diff = zeros(size(t));
    out.dC_nuc = zeros(size(t));
    out.mg = zeros(size(t));
    out.M1 = zeros(size(t));
    out.M2 = zeros(size(t));
    out.M3 = zeros(size(t));
    out.mu = Par.mu * ones(size(t));
    out.dPsat = zeros(size(t));
    out.rho_g = NaN(size(t));
    out.drdt = zeros(size(t));
    out.porosity = zeros(size(t));
    out.M = zeros(size(t));
    out.D_cin = NaN(size(t));

end


function [dydt,out]=before_fragmentation(t,y,Par)

global fw_interpolant psat_interpolant 

out.t = t;
out.mg = [y(1,:)]';

out.M0 = [y(2,:)]';
out.M1 = [y(3,:)]';
out.M2 = [y(4,:)]';
out.M3 = [y(5,:)]';
out.pm = [y(6,:)]';
out.Z  = [y(7,:)]';


% Conduit radius
out.a = Par.a * ones(size(t));      % Conduit radius
out.dadz = zeros(size(t));          % Conduit radius derivative in respect to z


out.Cm = Par.C0 - out.mg/Par.rho_melt;

%=============================================
% Nucleation rate as a function of supersaturation pressure

out.D = DH2Orhyolite(out.Cm,out.pm,Par.T);    %Diffusivity of volatile in melt [m^2/s]
% out.psat = find_psat(out.Cm,Par.T,Par.Xc);



out.psat = psat_interpolant(out.Cm);
out.pb = findPb(out.pm,out.psat,Par.T,fw_interpolant);  % Pressure in a nucleating bubble (Cluzel et al 2008)
out.sigma = SurfaceTension(out,Par.T); 

out.rc = 2*out.sigma ./ (out.pb-out.pm);
n0 = out.Cm * Par.rho_melt * Par.AV / Par.Mw;
a0 = n0.^(-1/3);
Vw = Molecular_VH2O(out.pm,Par.T);
A = (2*Vw .* n0.^2 .* out.D) ./ (a0) .* sqrt(out.sigma/(Par.KB*Par.T));
out.Wcl = (16*pi*out.sigma.^3) ./ (3*(out.pb-out.pm).^2);
out.J = A .* exp(-out.Wcl./(Par.KB*Par.T) * Par.ST_coeff);
%=============================================

% Equation of state for bubbles
r = out.M1 ./ out.M0;
% V = 4/3 * pi * r.^3;
V = 4/3 * pi * out.M3./out.M0;
m = out.mg./out.M0;
rho = m./V;
out.rho_g = rho;
[out.pg,  Kg] = EoS_H2O(rho,Par.T);


%=============================================
% Diffusion rate and bubble growth
out.cin = solubility(out.pg,Par.Xc,Par.T);
out.D_cin = DH2Orhyolite(out.cin,out.pg,Par.T);

mu = viscosity(Par.composition,Par.T,out.cin);
drdt = (r./(4*mu)) .* (out.pg-out.pm-2*out.sigma./r);
out.drdt = (drdt .* out.M0 + out.J .* out.rc)./out.M0;

mc = out.pb .* 4/3 * pi .* (out.rc).^3 / (Par.Rgas/Par.Mw * Par.T);



% out.dC_diff = 4*pi*r .* out.D_cin *Par.rho_melt .* (out.Cm-out.cin) .* out.M0; % dmg/dt
out.dC_nuc = mc .* out.J;
out.dC_diff = (4*pi*out.M2) .* out.D_cin *Par.rho_melt .* (out.Cm-out.cin)./r; % dmg/dt


dPsat = Psat_rate(out.psat,Par.T,(out.dC_diff + out.dC_nuc),Par);
out.dPsat = dPsat;

%======================================================================
% Magma decompression
out.mu = viscosity(Par.composition,Par.T,out.Cm);

out.porosity = out.M3 ./ (out.M3 + 3/(4*pi));
out.rho_magma = out.porosity .* out.rho_g +...
    (1-out.porosity) .* Par.rho_melt;
%----------------
K_magma = 1./((out.porosity./Kg) +...
    (1-(out.porosity)) ./ Par.K_melt);
% out.Kg = Kg;
% out.K_magma = K_magma;
C = (K_magma ./ out.rho_magma).^(1/2);
% C = max((K_magma ./ out.rho_magma).^(1/2), (Kg ./ out.rho_g).^(1/2));

out.C = C;
%---------------

U_ave = Par.Q./(out.rho_magma*pi.*out.a.^2);
out.U = U_ave;
M = out.U ./ C;
out.M = M;



% Friction factor
fric = 16*out.mu ./ (out.rho_magma.*U_ave*2.*out.a) + Par.f0;

% Decompression from Mastin et al 2000
    dPdz = (out.rho_magma * Par.g ...
     + out.rho_magma .* U_ave.^2 .* fric./out.a ...
     - 2 * out.rho_magma .* U_ave.^2 ./ out.a .* out.dadz)./(1-M.^2);
 
    out.dPdt = dPdz .* out.U;
    





%======================================================================
dmg  = out.dC_diff + out.dC_nuc;
dM0  = out.J;
dM1  = drdt .* out.M0 + out.J .* out.rc;
dM2  = 2 * drdt .* out.M1 + out.J .* (out.rc).^2;
dM3  = 3 * drdt .* out.M2 + out.J .* (out.rc).^3;
dPdt = -out.dPdt;
dzdt  = -out.U;



dydt = [dmg; dM0; dM1; dM2; dM3; dPdt; dzdt];

end


function [dydt,out]=after_fragmentation(t,y,Par)

global psat_interpolant 

out.t = t;
out.mg = [y(1,:)]';

out.M0 = [y(2,:)]';
out.M1 = [y(3,:)]';
out.M2 = [y(4,:)]';
out.M3 = [y(5,:)]';
out.pm = [y(6,:)]';
out.Z  = [y(7,:)]';
out.pg = [y(8,:)]';
out.a = [y(9,:)]';


% Conduit radius
out.a = Par.a * ones(size(t));      % Conduit radius
out.dadz = zeros(size(t));          % Conduit radius derivative in respect to z


out.Cm = Par.C0 - out.mg/Par.rho_melt;

%=============================================
out.D = DH2Orhyolite(out.Cm,out.pm,Par.T);    %Diffusivity of volatile in melt [m^2/s]
% out.psat = find_psat(out.Cm,Par.T,Par.Xc);

out.psat = psat_interpolant(out.Cm);
out.pb = out.psat;
out.sigma = NaN(size(t));

out.rc = NaN(size(t));



out.Wcl = NaN(size(t));
out.J = zeros(size(t));
%=============================================
% EoS within pyroclasts
vg_in = 4/3 * pi * out.M3;
[rho_in,Kg_in] = EoS_H2O_2(out.pg,Par.T);
mg_in = rho_in.*vg_in;


mg_out = out.mg - mg_in;
pg_out = out.pm;

% Gas properties
[rho_out,Kg_out] = EoS_H2O_2(pg_out,Par.T);
vg_out = mg_out ./ rho_out;



r = out.M1 ./ out.M0;
out.rho_g = mg_out./vg_out;

%=============================================
out.cin = solubility(out.pg,Par.Xc,Par.T);
out.D_cin = DH2Orhyolite(out.cin,out.pg,Par.T);

out.drdt = zeros(size(t));



% out.dC_diff = 4*pi*r .* out.D_cin *Par.rho_melt .* (out.Cm-out.cin) .* out.M0; % dmg/dt
out.dC_nuc = zeros(size(t));
out.dC_diff = (4*pi*out.M2) .* out.D_cin *Par.rho_melt .* (out.Cm-out.cin)./r; % dmg/dt


dPsat = Psat_rate(out.psat,Par.T,(out.dC_diff + out.dC_nuc),Par);
out.dPsat = dPsat;

%======================================================================
out.mu = viscosity(Par.composition,Par.T,out.Cm);

porosity_in = (vg_in) ./ (1+vg_in+vg_out);
porosity_out = (vg_out) ./ (1+vg_in+vg_out);
out.porosity = (vg_in+vg_out) ./ (1+vg_in+vg_out);
out.rho_magma = porosity_in .* rho_in ...
            + porosity_out .* rho_out ...
            + (1-porosity_in-porosity_out) * Par.rho_melt;
%----------------
K_magma = 1./((porosity_in./Kg_in) +...
    (porosity_out./Kg_out) + ...
    (1-porosity_in-porosity_out) ./ Par.K_melt);
C = (K_magma ./ out.rho_magma).^(1/2);
% C = max((K_magma ./ out.rho_magma).^(1/2), (Kg_out ./ rho_out).^(1/2));

out.C = C;
%---------------

U_ave = Par.Q./(out.rho_magma*pi.*out.a.^2);
out.U = U_ave;
M = U_ave ./ C;
Mc = 1;



fric = Par.f0;
% out.dadz(M>Mc) = 1/2 * (out.a(M>Mc)*Par.g./U_ave(M>Mc).^2 + fric);


    dPdz = (out.rho_magma * Par.g ...
     + out.rho_magma .* U_ave.^2 .* fric./out.a ...
     - 2 * out.rho_magma .* U_ave.^2 ./ out.a .* out.dadz)./(1-M.^2);
 
out.dPdt = dPdz .* out.U; 



out.M = M;



%======================================================================
dmg  = out.dC_diff + out.dC_nuc;
dM0  = zeros(size(t));
dM1  = zeros(size(t));
dM2  = zeros(size(t));
dM3  = zeros(size(t));
dPdt = -out.dPdt;
dzdt  = -out.U;


l1 = 1e-3;
permeability = 1e-13;
mu_gas = 1e-5;
Beta = 1./out.pg;
porosity_pyroclast = out.M3 ./ (3/4/pi + out.M3);
c = permeability./(mu_gas*porosity_pyroclast.*Beta);
tau_g = l1^2./c;
dpg  = (out.pm - out.pg)./tau_g;



dydt = [dmg; dM0; dM1; dM2; dM3; dPdt; dzdt; dpg; out.dadz];

end




function [dydt,out]=flaring(t,y,Par)

global psat_interpolant 

out.t = t;
out.mg = [y(1,:)]';

out.M0 = [y(2,:)]';
out.M1 = [y(3,:)]';
out.M2 = [y(4,:)]';
out.M3 = [y(5,:)]';
out.pm = [y(6,:)]';
out.Z  = [y(7,:)]';
out.pg = [y(8,:)]';
out.a = [y(9,:)]';





out.Cm = Par.C0 - out.mg/Par.rho_melt;

%=============================================
out.D = DH2Orhyolite(out.Cm,out.pm,Par.T);    %Diffusivity of volatile in melt [m^2/s]
% out.psat = find_psat(out.Cm,Par.T,Par.Xc);

out.psat = psat_interpolant(out.Cm);
out.pb = out.psat;
out.sigma = NaN(size(t));

out.rc = NaN(size(t));



out.Wcl = NaN(size(t));
out.J = zeros(size(t));
%=============================================
% EoS within pyroclasts
vg_in = 4/3 * pi * out.M3;
[rho_in,Kg_in] = EoS_H2O_2(out.pg,Par.T);
mg_in = rho_in.*vg_in;


mg_out = out.mg - mg_in;
pg_out = out.pm;

% Gas properties
[rho_out,Kg_out] = EoS_H2O_2(pg_out,Par.T);
vg_out = mg_out ./ rho_out;



r = out.M1 ./ out.M0;
out.rho_g = mg_out./vg_out;

%=============================================
out.cin = solubility(out.pg,Par.Xc,Par.T);
out.D_cin = DH2Orhyolite(out.cin,out.pg,Par.T);

out.drdt = zeros(size(t));



% out.dC_diff = 4*pi*r .* out.D_cin *Par.rho_melt .* (out.Cm-out.cin) .* out.M0; % dmg/dt
out.dC_nuc = zeros(size(t));
out.dC_diff = (4*pi*out.M2) .* out.D_cin *Par.rho_melt .* (out.Cm-out.cin)./r; % dmg/dt


dPsat = Psat_rate(out.psat,Par.T,(out.dC_diff + out.dC_nuc),Par);
out.dPsat = dPsat;

%======================================================================
out.mu = viscosity(Par.composition,Par.T,out.Cm);

porosity_in = (vg_in) ./ (1+vg_in+vg_out);
porosity_out = (vg_out) ./ (1+vg_in+vg_out);
out.porosity = (vg_in+vg_out) ./ (1+vg_in+vg_out);
out.rho_magma = porosity_in .* rho_in ...
            + porosity_out .* rho_out ...
            + (1-porosity_in-porosity_out) * Par.rho_melt;
%----------------
K_magma = 1./((porosity_in./Kg_in) +...
    (porosity_out./Kg_out) + ...
    (1-porosity_in-porosity_out) ./ Par.K_melt);
C = (K_magma ./ out.rho_magma).^(1/2);
% C = max((K_magma ./ out.rho_magma).^(1/2), (Kg_out ./ rho_out).^(1/2));

out.C = C;
%---------------

U_ave = Par.Q./(out.rho_magma*pi.*out.a.^2);
out.U = U_ave;
M = U_ave ./ C;



fric = Par.f0;

out.dadz = 1/2 * (out.a*Par.g./C.^2 + fric);
dadt = out.dadz .* out.U;

 


dPdz = (out.rho_magma * Par.g ...
     + out.rho_magma .* U_ave.^2 .* fric./out.a ...
     - 2 * out.rho_magma .* U_ave.^2 ./ out.a .* out.dadz)./(1-M.^2);
 
out.dPdt = dPdz .* out.U; 



out.M = M;



%======================================================================
dmg  = out.dC_diff + out.dC_nuc;
dM0  = zeros(size(t));
dM1  = zeros(size(t));
dM2  = zeros(size(t));
dM3  = zeros(size(t));
dPdt = -out.dPdt;
dzdt  = -out.U;


l1 = 10e-2;
permeability = 1e-13;
mu_gas = 1e-5;
Beta = 1./out.pg;
porosity_pyroclast = out.M3 ./ (3/4/pi + out.M3);
c = permeability./(mu_gas*porosity_pyroclast.*Beta);
tau_g = l1^2./c;
dpg  = (out.pm - out.pg)./tau_g;



dydt = [dmg; dM0; dM1; dM2; dM3; dPdt; dzdt; dpg; dadt];

end







function dPsat = Psat_rate(Pw,T,dmg,Par)
% The rate of decrease in water saturation pressure
    a1 = 354.94;
    a2 = 9.623;
    a3 = -1.5223;
    a4 = 0.0012439;
    a5 = -1.084e-4;
    a6 = -1.362e-5;
    
    
    
    Pw = Pw/1e6;
    
    coeff = (1/2*a1*Pw.^(-1/2) + a2 + 3/2*a3*Pw.^(1/2))/T + 3/2*a4*Pw.^(1/2);
    
    dC = dmg/Par.rho_melt*100;
    dPsat = dC./coeff * 1e6;

end


function [Cw,Cc] = solubility(P,Xc,T)
% Liu et al 2005
    Xw = 1 - Xc;
    Pw = Xw .* P/1e6;
    Pc = Xc .* P/1e6;
  
    a1 = 354.94;
    a2 = 9.623;
    a3 = -1.5223;
    a4 = 0.0012439;
    a5 = -1.084e-4;
    a6 = -1.362e-5;
    %
    b1 = 5668;
    b2 = 0.4133;
    b3 = 2.041e-3;
    b4 = -55.99;
    
    
    Cw = (a1*Pw.^(1/2) + a2*Pw + a3*Pw.^(3/2))./T + a4*Pw.^(3/2) +...
      Pc.*(a5*Pw.^(1/2) + a6*Pw); %H2O content in wt.%
    Cc = b1*Pc/T + Pc.*(b2*Pw.^(1/2) + b3*Pw.^(3/2)) + b4*Pc.*Pw/T;
  
    Cw = Cw/100; % fraction
    Cc = Cc/1e6; % fraction
    
end


function psat = find_psat(cm,T,Xc)
% Inverse of solubility, saturation pressure as a function of water
% concentration
global Pinitial

psat = zeros(size(cm));


p2 = Pinitial;
p1 = 1e5;


% if cm < .06
%     keyboard
% end

for i = 1:length(cm)
 
  fun =  @(x) abs(solubility(x,Xc,T) - cm(i)); 
  
 psat(i) = fminbnd(fun,p1,p2,optimset('TolX',1e-2));   
    
    
end


end


function eta = viscosity(composition,T,H2O)

% Viscosity, Hui & Zhang 2007

% SiO2 = 76.53e-2;
% TiO2 = .06e-2;
% Al2O3 = 13.01e-2;
% FeO = .79e-2;
% MnO = .08e-2;
% MgO = .02e-2;
% CaO = .74e-2;
% Na2O = 3.87e-2;
% K2O = 4.91e-2;

mO      = 15.9994;
mSiO2   = 28.0855 + 2*mO;
mTiO2   = 47.88 + 2*mO;
mAl2O3  = 2*26.98154 + 3*mO;
mFeO    = 55.847 + mO;
mFe2O3  = 2*55.847 + 3*mO; 
mMnO    = 54.9380 + mO;
mMgO    = 24.305 + mO;
mCaO    = 40.08 + mO;
mNa2O   = 2*22.98977 + mO;
mK2O    = 2*39.0983 + mO;
mH2O    = 2*1.00794 + mO;

nSiO2 = composition.SiO2 / mSiO2;
    nTiO2 = composition.TiO2 / mTiO2;
    nAl2O3 = composition.Al2O3 / mAl2O3;
    nFeO = composition.FeO / mFeO;
    nFe2O3 = composition.Fe2O3 / mFe2O3; % Added by CR, Sep 2023
    nMnO = composition.MnO / mMnO;
    nMgO = composition.MgO / mMgO;
    nCaO = composition.CaO / mCaO;
    nNa2O = composition.Na2O / mNa2O;
    nK2O = composition.K2O / mK2O;
%    nP2O5 = P2O5 / mP2O5;
    nH2O = H2O / mH2O;
    nFeMnO = nFeO + nMnO + nFe2O3; % Fe2O3 included by CR, Sep 2023, assuming HZ(2007) FeO_tot wt.% is equivalent to converted molar mass...
    
    nNaK = 2*nNa2O + 2*nK2O;
    nAl = 2*nAl2O3;
    if nNaK <= nAl
        nNaKAlO2 = nNaK;
        nAl2O3ex = (nAl-nNaK)/2;
        nNaK2Oex = 0;
    else
        nNaKAlO2 = nAl;
        nNaK2Oex = (nAl-nNaK)/2;
        nAl2O3ex = 0;
    end
    
    nmol = nSiO2 + nTiO2 + nFeMnO + ...
        nMgO + nCaO + nNaKAlO2 + nAl2O3ex + nNaK2Oex + nH2O;
    %
    XSiO2 = nSiO2./nmol;
    XTiO2 = nTiO2./nmol;
    XFeMnO = nFeMnO./nmol;
    XMgO = nMgO./nmol;
    XCaO = nCaO./nmol;
 %   XP2O5 = nP2O5./nmol;
 XP2O5 = 0;
    XNaKAlO2 = nNaKAlO2./nmol;
    XAl2O3ex = nAl2O3ex./nmol;
    XNaK2Oex = nNaK2Oex./nmol;
    XH2O = nH2O./nmol;
    %
    Z = XH2O.^(1./(1 + 185.797./T));
    
    eta = 10.^( ...
        ( ...
        -  6.83*XSiO2 ...
        - 170.79*XTiO2 ...
        - 14.71*XAl2O3ex ...
        - 18.01*XMgO ...
        - 19.76*XCaO ...
        - 34.31*XNaK2Oex ...
        - 140.38*Z ...
        + 159.26*XH2O ...
        - 8.43*XNaKAlO2 ...
        ) ...
        + ...
        ( ...
        + 18.14*XSiO2 ...
        + 248.93*XTiO2 ...
        + 32.61*XAl2O3ex ...
        + 25.96*XMgO ...
        + 22.64*XCaO ...
        - 68.29*XNaK2Oex ...
        + 38.84*Z ...
        - 48.55*XH2O ...
        + 16.12*XNaKAlO2 ...
        )* 1000./T ...
        + ...
        exp( ...
        ( ...
        + 21.73*XAl2O3ex ...
        - 61.98*XFeMnO ...
        - 105.53*XMgO ...
        - 69.92*XCaO ...
        - 85.67*XNaK2Oex ...
        + 332.01*Z ...
        - 432.22*XH2O ...
        - 3.16*XNaKAlO2 ...
        ) ...
        + ...
        ( ...
        + 2.16*XSiO2 ...
        - 143.05*XTiO2 ...
        - 22.10*XAl2O3ex ...
        + 38.56*XFeMnO ...
        + 110.83*XMgO ...
        + 67.12*XCaO ...
        + 58.01*XNaK2Oex ...
        + 384.77*XP2O5...
        - 404.97*Z ...
        + 513.75*XH2O ...
        )*1000./T ...
        ) ...
        );
end


function D = DH2Orhyolite(concentration,pm,T)
% Diffusivity
diffusivity = 'Zhang_Behrens_2000';
%
if strcmp(diffusivity,'Zhang_Behrens_2000')
    %H2O bulk diffusivity (D in m^2/s) in rhyolitic melt after Zhang and Behrens (2000)
    %
    %T is in Kelvin
    P = pm/1e6; %from Pa to MPa
    C = 100*concentration; % from fractional to percent
    X = (C/18.015) ./ (C/18.015 + (100-C)/32.49); %mole fraction of total H2O on a single oxygen basis
    m = -20.79 - 5030./T - 1.4*P./T;

    D = 1e-12 .* X .* exp(m) .* ( ...
        1 +exp( ...
        56 + m + X.*(-34.1 + 44620./T + 57.3*P./T) ...
        - sqrt(X).*(0.091 + 4.77e6./T.^2) ...
        ) ...
        );
elseif strcmp(diffusivity,'Ni_Zhang_2008')
    %H2O bulk diffusivity (D in m^2/s) in rhyolitic melt after
    %Ni and Zhang, Chemical Geology, 250, 68-78 (2008)
    %
    %T is in Kelvin
    P = pm/1e9; %from Pa to GPa
    C = 100*concentration; % from fractional to percent
    X = (C/18.015) ./ (C/18.015 + (100-C)/32.49); %mole fraction of total H2O on a single oxygen basis

    D = 1e-12 * X.*exp( ...
        13.47 - 49.996*X + 7.0827*sqrt(X) + 1.8875*P ...
        - (9532.3 - 91933*X + 13403*sqrt(X) + 3625.6*P)./T ...
        );


end



end


function pb = findPb(pm,psat,T,fw_interpolant)
% Pressure in a bubble nucleus (Cluzel et al 2008)

KB = 1.38e-23;                      %Boltzman constant [J/K]

pb = zeros(size(pm));


for i = 1:length(pm)
    Vw = Molecular_VH2O(pm(i),T);
const = exp(Vw/KB/T*(pm(i)-psat(i)))...
    .* fw_interpolant(psat(i));

p1 = pm(i);
p2 = psat(i);

fun =  @(x) abs(fw_interpolant(x) - const); 
pb(i) = fminbnd(fun,p1,p2,optimset('TolX',1e-1)); 


end



end


function [Pg, Kg] = EoS_H2O(rho,T)

% Modified Redlich and Kwong EoS for water vapor from Holloway 1977


R = 83.12;          % Gas constant cm^3.bar/(deg mole)
M = 18.01528e-3;    % Molar mass of water kg/mol

TC = T - 273.15;    % Degree Cel

ao = 35e6;
b = 14.6;
a = 166.8e6 - 193080*TC + 186.4*TC.^2 - 0.071288*TC.^3;



V = (1./rho)*M*1e6;      % molar volume (cm^3/mol)

Pg = R*T./(V-b) - a./(sqrt(T).*V.*(V+b));

Pg = Pg * 1e5;          % bar to Pa


% Bulk Modulus (V is in molar volume)
dpdv = -R*T./(V-b).^2 + a./(sqrt(T)) * (2*V+b)./(V.*(V+b)).^2;
dvdrho = -(1./rho).^2 * M * 1e6;

k = 1; %isothermal (1.33 for isenthropic)
Kg = k * rho .* dpdv .* dvdrho; % 1.33 is Cp/Cv
Kg = Kg * 1e5;

Kg(Pg==0) = inf;

end


function [rho, Kg] = EoS_H2O_2(Pg,T)

% Modified Redlich and Kwong EoS for water vapor from Holloway 1977


R = 83.12;          % Gas constant cm^3.bar/(deg mole)
M = 18.01528e-3;    % Molar mass of water kg/mol

TC = T - 273.15;    % Degree Cel

ao = 35e6;
b = 14.6;
a = 166.8e6 - 193080*TC + 186.4*TC.^2 - 0.071288*TC.^3;

Pg = Pg / 1e5;          % Pa to bar

rho = NaN(size(Pg));

for i = 1:length(Pg)
    p = [Pg(i) -R*T a./sqrt(T)-Pg(i)*b^2-R*T*b -a./sqrt(T)*b];
    for j = 1:4
        if isnan(p(j)) || isinf(p(j))
            keyboard
        end
    end
    r = roots(p);
    V = r(imag(r)==0&real(r)>0);
    rho(i) = (1/V)*M*1e6;      % molar volume (cm^3/mol)
end

V = (1./rho)*M*1e6; 

% Bulk Modulus
dpdv = -R*T./(V-b).^2 + a./(sqrt(T)) * (2*V+b)./(V.*(V+b)).^2;
dvdrho = -(1./rho).^2 * M * 1e6;

k = 1; %isothermal (1.33 for isenthropic)
Kg = k * rho .* dpdv .* dvdrho; % 1.33 is Cp/Cv
Kg = Kg * 1e5;

Kg(Pg==0) = inf;

end


function Vw = Molecular_VH2O(P,T)
% Molar volume of water (Ochs & Lange, Science 1999)

AV = 6.022e23;      % Avogadro number
T = T - 273.15;     % K to C
P =  P/1e5;         % Pa to bar


V = 22.9 + 9.5e-3*(T-1000) - 3.2e-4*(P-1);  % molar volume of water (cm^3)

Vw = V * 1e-6 / AV;

end

function ST = SurfaceTension(Par,T)
% Hajimirza et al. JGR 2019

delta = .3200e-3;   % micro meter
alpha = .51;

Psat = Par.psat/1e6;    % Pa to MPa
Pb = Par.pb/1e6;
Pm = Par.pm/1e6;    % Pa to MPa
T = T - 273.15;         % K to C

ST_B = 1.21e-1 * exp(-2.24e-2*Psat)...
     + 1.47e-1 * exp(-1.90e-3*Psat)...
     + 7.5e-5 * (T-1000);
 
ST_inf = ST_B * (1-alpha);

ST = ST_inf + delta * (Pb-Pm);

end







function [position,isterminal,direction] = EventsFcn1(t,y,Par)
% Event function for before nucleation: ODE solver stops when either of the following happens:
% Pressure = 1 atm
% Bubble number density = Minimum threshold for bubble number density

position = [log10(y(1))-log10(Par.ND) y(3)]; % The value that we want to be zero
isterminal = [1 1];  % Halt integration 
direction = [1 -1];   % Event function increasing
end

function [position,isterminal,direction] = EventsFcn2(t,y,Par)
% Event function for before fragmentation: ODE solver stops when either of the following happens:
% Pressure = 1 atm
% Z = 0
% Porosity = Fragmentation porosity


phi = y(5)/(3/4/pi + y(5));

position = [Par.phi_frag-phi y(6)-Par.pf y(7)]; % The value that we want to be zero
isterminal = [1 1 1];  % Halt integration 
direction = [-1 -1 -1];   % Event function increasing
end

function [position,isterminal,direction] = EventsFcn3(t,y,Par)
% Event function for after fragmentation: ODE solver stops when either of the following happens:
% Pressure = 1 atm
% Z = 0
% Mach number = 0.99

mg = y(1);
M3 = y(5);
pm = y(6);
pg = y(8);
Cm = Par.C0 - mg/Par.rho_melt;



vg_in = 4/3 * pi * M3;
[rho_in,Kg_in] = EoS_H2O_2(pg,Par.T);
mg_in = rho_in.*vg_in;


mg_out = mg - mg_in;
pg_out = pm;

% Gas properties
[rho_out,Kg_out] = EoS_H2O_2(pg_out,Par.T);
vg_out = mg_out ./ rho_out;



% rho_g = mg_out./vg_out;

porosity_in = (vg_in) ./ (1+vg_in+vg_out);
porosity_out = (vg_out) ./ (1+vg_in+vg_out);
rho_magma = porosity_in .* rho_in ...
            + porosity_out .* rho_out ...
            + (1-porosity_in-porosity_out) * Par.rho_melt;
K_magma = 1./((porosity_in./Kg_in) +...
    (porosity_out./Kg_out) + ...
    (1-porosity_in-porosity_out) ./ Par.K_melt);
% C = max((K_magma ./ rho_magma).^(1/2), (Kg_out ./ rho_out).^(1/2));

C = (K_magma ./ rho_magma).^(1/2);


U_ave = Par.Q./(rho_magma*pi.*Par.a^2);
M = U_ave ./ C;

position = [y(6)-Par.pf y(7) Par.Mc-M]; % The value that we want to be zero
isterminal = [1 1 1];  % Halt integration 
direction = [-1 -1 -1]; % The zero is approached from up to down
end


function [position,isterminal,direction] = EventsFcn4(t,y,Par)
% Event function for after fragmentation: ODE solver stops when either of the following happens:
% Pressure = 1 atm
% Z = 0
% Mach number = 0.99

mg = y(1);
M3 = y(5);
pm = y(6);
pg = y(8);
a = y(9);
Cm = Par.C0 - mg/Par.rho_melt;



vg_in = 4/3 * pi * M3;
[rho_in,Kg_in] = EoS_H2O_2(pg,Par.T);
mg_in = rho_in.*vg_in;


mg_out = mg - mg_in;
pg_out = pm;

% Gas properties
[rho_out,Kg_out] = EoS_H2O_2(pg_out,Par.T);
vg_out = mg_out ./ rho_out;



% rho_g = mg_out./vg_out;

porosity_in = (vg_in) ./ (1+vg_in+vg_out);
porosity_out = (vg_out) ./ (1+vg_in+vg_out);
rho_magma = porosity_in .* rho_in ...
            + porosity_out .* rho_out ...
            + (1-porosity_in-porosity_out) * Par.rho_melt;
K_magma = 1./((porosity_in./Kg_in) +...
    (porosity_out./Kg_out) + ...
    (1-porosity_in-porosity_out) ./ Par.K_melt);
% C = max((K_magma ./ rho_magma).^(1/2), (Kg_out ./ rho_out).^(1/2));

C = (K_magma ./ rho_magma).^(1/2);


U_ave = Par.Q./(rho_magma*pi.*a^2);
M = U_ave ./ C;

position = [y(6)-Par.pf y(7) .999-M]; % The value that we want to be zero
isterminal = [1 1 1];  % Halt integration 
direction = [-1 -1 -1]; % The zero is approached from up to down
end




