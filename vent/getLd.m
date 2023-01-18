function [D,pO] = getLd(C,P,model)
% Return jet properties after vent decompression, including decompression
% vertical length scale.
%
% C     = struct of vent params. Required fields:
%           rho_magma   mixture density
%           rho_g       gas density
%           rho_s       solids bulk density
%           C           mixture sound speed
%           u           mixture velocity
%           n_0         free gas mass fraction
%           a           vent radius
%           pm          vent pressure
%           pf          vent ambient pressure
%           T           mixture temperature
%           Ks          total bulk modulus of solids
%           Zw          water depth above vent (m)
%
% P      = plumeSource struct (see getPlumeSource)
%
% model = decompression model
%           (1) Free decompression model following Woods and Bower 1995, and
%               Colucci et al 2014. Decompression pressure is estimated
%               from ambient plus dynamic pressure as P_d + rho*c^2/2 =
%               P_a. L_d ~ 2*(r_d - r_v)
%           (2) Same as (1), but with decompression velocity limited to
%               M<=1, similar to Ogden 2008
%           (3) Same as (1), but with pressure averaged between vent
%               pressure and pressure at L_d
%           (4) Same as (2), but with pressure averaged between vent
%               pressure and pressure at L_d
%           (5) Ogden 2008 Mach disk height scaling, with M=1 above disk
%           (6) Same as (5), but with pressure averaged between vent
%               pressure and pressure at L_d
%
% RETURNS: D = struct of decompression parameters with fields:
%               Ld      Decompression length        (m above vent)
%               pd      Internal jet pressure at Ld (Pa)
%               pf      Ambient pressure at Ld      (Pa)
%               u       Jet velocity at Ld          (m/s)
%               c       Gas sounds speed at Ld      (m/s)
%               rho_d   Jet bulk density at Ld      (kg/m^3)
%               rho_g   Jet gas density at Ld       (kg/m^3)
%               Kg      Gas bulk modulus at Ld      (Pa)
%               chi      Gas volume fraction at Ld   --
%               a       Jet radius at Ld            (m)
%               
%           pO = updated plume source struct
% 
% C Rowell, May 2021
%

if nargin<3
    model = 2; % Default decompression model
end

pO = P;
% Constants
P.gamma     = heatCapacityP(C.rho_g,C.T)/heatCapacityV(C.rho_g,C.T); %1.33; % adiabatic index for water
P.g         = 9.81; % g
P.LdScale   = 1;    % Blast Ld constant
P.Rw        = 461;  % Vapour gas constant
P.rho_l     = 1000; % density water

% P.c_air     = 343;  % sound speed air
% P.rho_air   = 1;    % density air
% P.c_water   = 1481; % sound speed water

% Get parameters
P.pa       = C.pf; % Ambient pressure at vent
P.Patmo    = C.pf - P.rho_l*9.81*C.Zw; % Pressure at air surface
% C.Zw       = Zw;    % Keep water depth
% P.P_d       = P.pa/(1+P.gamma/2); % Initial guess

%%

    if ~P.useDecompressLength
        D = Ldmodel(C,P,9);
        
    elseif C.Zw ==0
        D = Ldmodel(C,P,7);
        
    else
        % Run first decompression
        D = Ldmodel(C,P,model);
        
        % May need to re-run here if D.pd was much greater than Patmo?
%         if D.Ld_gt_Zw && D.pd~=D.pf
%             D2 = getLd(C,P,7);
%         end
    end
        
        % Check Ld vs Ldmin? Not needed for now, Ld is bounded at current Ldmin
       
        % Readjust for pressure at Ld
        % --> Let bubble gas equilibrate for now
%         n_bs0 = sum(P.nsi.*P.ni); % Initial mass frac of bubbles in pyroclasts
%         [D.nsi,D.ni,D.Rgsd,D.rhoi,D.pori] = getPSD(P.rho_m, P.phi0, D.pf, C.T, P.D, P.phiSz_min, P.phiSz_max);
%         n_bs  = sum(D.nsi.*D.ni);
%         C.n_0 = ((1-P.n_0)*(1-n_bs0) - 1 + n_bs)./(n_bs - 1);
% %         D.n_0 = C.n_0;
%         
%         % --> Currently NOT accounting for dynamic pressure
%         D = paramsAtPd(D.pf,C,D);
%         D.u = D.c;
%         D.a = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2);

        % Check Ld vs Zw 
%         if D.Ld_gt_Zw
%             % Get new params and run second decomp
%         end
%     end
    
    D.model      = model;
    % Copy over relevant plume params
    pO.Ld        = D.Ld;
    pO.vh0       = pO.vh0 + pO.Ld;
    pO.n_0       = D.n_0;
    pO.ni        = D.ni;
    pO.nsi       = D.nsi;
    pO.rhoi      = D.rhoi;
    pO.Pg        = D.pf;
    pO.r_0       = D.a; %r_out;
    pO.u_0       = D.u;
    pO.rho_B0    = D.rho_d; %rhoB_out;
    pO.rho_g0    = D.rho_g; %rho_gd;
    
     %% === Temp plots to run some tests
    if 0==1
        %         searchFun = @(Pd) LdSearchFun(Pd,C,P,model);
%         Pf = fminsearch(searchFun,P.P_d);
%         P.P_d = Pf;
%         D = Ldmodel(C,P,model);
%         Ld = D.Ld; 
        %         freeLdtests
        
        ZwT = linspace(C.Zw,0,101);
        Pz = linspace(P.pa,P.Patmo,101);
        Z  = linspace(0,C.Zw,101);
            
        Ld_base = zeros(length(Z),3);
        Ld_iter = Ld_base;
        a_base  = Ld_base;
        a_iter  = Ld_base;
        u_base  = Ld_base;
        u_iter  = Ld_base;
        c_base  = Ld_base;
        c_iter  = Ld_base;

%         D0 = Ldmodel(C,P,model);
        P0 = P;
        pf1 = zeros(size(Z));
%         D2 = D0;
        P1 = P0;
%         P2 = P0;
%         P3 = P0;
%         P4 = P0;
        for ii = 1:length(Pz)
            P1.P_a = Pz(ii);
%             P1.P_d = Pz(ii);
            P1.Zw   = ZwT(ii);
            
            D1 = Ldmodel(C,P1,1);
            Ld_base(ii,1) = D1.Ld;
            a_base(ii,1)  = D1.a;
            u_base(ii,1)  = D1.u;
            c_base(ii,1)  = D1.c;
            pf1(ii) = D1.P_f;
            
            D2 = Ldmodel(C,P1,2);
            Ld_base(ii,2) = D2.Ld;
            a_base(ii,2)  = D2.a;
            u_base(ii,2)  = D2.u;
            c_base(ii,2)  = D2.c;
            
            D3 = Ldmodel(C,P1,5);
            Ld_base(ii,3) = D3.Ld;
            a_base(ii,3)  = D3.a;
            u_base(ii,3)  = D3.u;
            c_base(ii,3)  = D3.c;
            
            D4 = Ldmodel(C,P1,3);
            Ld_iter(ii,1) = D4.Ld;
            a_iter(ii,1)  = D4.a;
            u_iter(ii,1)  = D4.u;
            c_iter(ii,1)  = D4.c;
             
            D5 = Ldmodel(C,P1,4);
            Ld_iter(ii,2) = D5.Ld;
            a_iter(ii,2)  = D5.a;
            u_iter(ii,2)  = D5.u;
            c_iter(ii,2)  = D5.c;
            
            D6 = Ldmodel(C,P1,6);
            Ld_iter(ii,3) = D6.Ld;
            a_iter(ii,3)  = D6.a;
            u_iter(ii,3)  = D6.u;
            c_iter(ii,3)  = D6.c;
            
%             Ud1(ii) = D1.u;
%             Pf1(ii) = max([(P1.P_a-D1.Ld*P.g*P.rho_l) P.Patmo]);
            
%             rho_z = (C.n_0*P.Rw*C.T./Pz(ii) + (1-C.n_0)./C.rho_s).^(-1); % Bulk density equilibrated to local pressure
%             cz    = sqrt(P.gamma.*Pz(ii)./rho_z); % Approx sound speed at balanced pressure
%             P2.P_d = Pz(ii)+(0.5).*cz.^2.*rho_z;
%             P2.P_a = Pz(ii);
%             P2.P_d = max([Pz(ii) P.Patmo])./(1+P.gamma/2); %max([Pz(ii)./(1-P.gamma/2) P.Patmo]);
%             P2.Zw   = ZwT(ii);
% %             rho_z = (C.n_0*P.Rw*C.T./P2.P_d + (1-C.n_0)./C.rho_s).^(-1); % Bulk density equilibrated to local pressure
% %             cz    = sqrt(P.gamma.*P2.P_d./rho_z); % Approx sound speed at balanced pressure
%             D2 = Ldmodel(C,P2,model);
%             Ld2(ii) = D2.Ld;
%             Ud2(ii) = D2.u;
%             Pf2(ii) = max([P2.P_a-Ld2(ii)*P.g*P.rho_l P.Patmo]);
         
            
            % Same curves but with 1 iteration on Pz
            
%             P3.Zw   = ZwT(ii);
%             P3.P_d = mean([Pz(ii) Pf1(ii)]);
%             D3 = Ldmodel(C,P3,model);
%             Ld3(ii) = D3.Ld;
%             Ud3(ii) = D3.u;
% %             Pf3(ii) = max([(P3.P_a-D3.Ld*P.g*P.rho_l) P.Patmo]);
%             
%             P4.Zw   = ZwT(ii);
% %             P4.P_d = Pf2(ii);
% %             rho_z = (C.n_0*P.Rw*C.T./P4.P_d + (1-C.n_0)./C.rho_s).^(-1); % Bulk density equilibrated to local pressure
% %             cz    = sqrt(P.gamma.*P4.P_d./rho_z); % Approx sound speed at balanced pressure
%             P4.P_d = max([mean([Pz(ii) Pf2(ii)]) P.Patmo])./(1+P.gamma/2);
%             D4 = Ldmodel(C,P4,model);
%             Ld4(ii) = D4.Ld;
%             Ud4(ii) = D4.u;          
%             Pf4(ii) = max([(P4.P_a-D4.Ld*P.g*P.rho_l) P.Patmo]);
           
        end
        ll{1} = '2dr, M<~1.8';
        ll{2} = '2dr, M<=1';
        ll{3} = 'Ogden';
        ll{4} = '2dr, M<~1.8, Piter';
        ll{5} = '2dr, M<=1, Piter';
        ll{6} = 'Ogden, Piter';
        ll{7} = 'pgZ, Zw=100';
        ll{8} = 'pgH_w';        

        figure
        plot(Pz/1e6,Ld_base,'LineWidth',2)
        hold on
        set(gca,'ColorOrderIndex',1)
        plot(Pz/1e6,Ld_iter,'--','LineWidth',2)
%         plot(Pf1/1e6,Ld1)
%         plot(Pf2/1e6,Ld2)
%         set(gca,'ColorOrderIndex',1)
%         plot(Pz/1e6,Ld3,'--','LineWidth',2)
%         plot(Pz/1e6,Ld4,'--','LineWidth',2)
        plot(Pz/1e6,Z,'--k','LineWidth',2)
        plot(Pz/1e6,ZwT,'k')
        plot(pf1/1e6,Ld_base(:,1),'g','LineWidth',2)
%         set(gca,'YDir','reverse')
        xlabel('P_a (MPa)')
        ylabel('Height (m a.v.l.)')
        grid on
        axis tight
        legend(ll) %{'$P_d = P_a$','$P_d = P_a + \rho u^2/2$','Hydrostatic'},'Interpreter','Latex')
%         title('L_d = 2(r_d - r_v)')

        % Check radii
        figure('position',[50 50 1000 600])
        subplot(1,2,1)
        plot(Pz/1e5,a_base,'LineWidth',2)
        hold on
        set(gca,'ColorOrderIndex',1)
        plot(Pz/1e5,a_iter,'--','LineWidth',2)
        xlabel('P_a (MPa)')
        ylabel('r_d (m)')
        grid on
        axis tight
        legend(ll{1:6})
        
        subplot(1,2,2)
        plot(Pz/1e5,u_base,'LineWidth',2)
        hold on
        plot(Pz/1e5,u_iter,'LineWidth',2)
        set(gca,'ColorOrderIndex',1)
        plot(Pz/1e5,c_base,'--','LineWidth',2)
        plot(Pz/1e5,c_iter,'--','LineWidth',2)
        xlabel('P_a (MPa)')
        ylabel('Velocity (m/s)')
        grid on
        axis tight
        legend([ll{1:6} {'Sound speed'}])

        % === End bullshit
    end

end


function [D] = Ldmodel(C,P,model)
% D struct fields: rho_B, rho_g, u, a, Ld, P_d

    switch model
        case 1 % Free decomp with gamma scaling up to M~1.75
%             D = freeDecomp(C,P);
            cScale       = 1+1./P.gamma-1./(P.gamma).*(P.pa./C.pm);
            D.pd        = P.pa./(1+P.gamma.*cScale.^2/2);                        % Initial estimate for final pressure
            D            = paramsAtPdFull(D.pd,C,D,P);
            D.u          = C.u + 1./(C.rho_magma.*C.u).*(C.pm-D.pd);   % Decompression velocity
            D.a          = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); % new radius
            D.Ld         = P.LdScale*2.*(D.a-C.a);                         % Decompression length scale
            
        case 2 % Free decomp, limited to M<=1
            D.pd        = P.pa./(1+P.gamma./2);                         % Initial estimate for final pressure using dynamic pressure scale
            D            = paramsAtPdFull(D.pd,C,D,P);
%             if (C.pm/C.pf)<1.05
%                 D.u = C.u + 1./(C.rho_magma.*C.u).*(C.pm-C.pf);
%             else
%                 D.u          = D.c;                                         % Decompression velocity
%             end
            D.u          = min([D.c (C.u + 1./(C.rho_magma.*C.u).*(C.pm-C.pf))]);
            D.a          = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); % new radius
            D.Ld         = P.LdScale*2.*(D.a-C.a);                         % Decompression length scale
            
        case 3 % One iteration on model 1 with averaged pressure
            D            = Ldmodel(C,P,1);
            cScale       = 1+1./P.gamma-1./(P.gamma).*(P.pa./C.pm);
            gScale       = (1+P.gamma.*cScale.^2/2);
            D.pd         = mean([D.pd (D.pf)./gScale]);
            D            = paramsAtPdFull(D.pd,C,D,P);
            D.u          = C.u + 1./(C.rho_magma.*C.u).*(C.pm-D.pd);   % Decompression velocity
            D.a          = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); % new radius
            D.Ld         = P.LdScale*2.*(D.a-C.a);                         % Decompression length scale            
            
        case 4 % One iteration on model 2 with averaged pressure
            D            = Ldmodel(C,P,2);
            gScale       = (1+P.gamma./2);
            D.pd         = mean([D.pd (D.pf)./gScale]); % Initial estimate for final pressure
            D            = paramsAtPdFull(D.pd,C,D,P);
            D.u          = D.c;                                            % Decompression velocity
            D.a          = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); % new radius
            D.Ld         = P.LdScale*2.*(D.a-C.a);                         % Decompression length scale
            
        case 5 % Ogden 2008 Mach disk height scale
%             D.pd        = P.pa./(1+P.gamma./2);                           % Initial estimate for final pressure
            D.pd         = P.pa;
            D            = paramsAtPdFull(D.pd,C,D,P);
            D.a          = C.a.*(C.pm./D.pd).^(1/2);                      % Ogden 2008 disk radius
            D.u          = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.a.^2) ); % new velocity
            D.Ld         = (pi*C.a.^2.*C.pm./D.pd).^(1/2);             % Ogden 2008 Mach disk height
            
        case 6 % Ogden 2008 Mach disk height scale with 1 iteration on averaged pressure
            D            = Ldmodel(C,P,5);
            D.pd         = mean([P.pa D.pf]);
            D            = paramsAtPdFull(D.pd,C,D,P);
            D.u          = D.c;                                         % Decompression velocity
            D.a          = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); % new radius
            D.Ld         = (pi*C.a.^2.*C.pm./D.pd).^(1/2);             % Ogden 2008 Mach disk height
        
        case 7 % Use decompression length, simple free decompression from Woods 1995/Colucci 2014
               % - (Equivalent to model 2 with no dynamic pressure adjustment). Decompression velocity limited to Mach 1
            D.pd        = P.pa;
            D           = paramsAtPdFull(D.pd,C,D,P);
%             D.u         = D.c;
            D.u          = min([D.c (C.u + 1./(C.rho_magma.*C.u).*(C.pm-C.pf))]);
            D.a         = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); % new radius
            D.Ld        = P.LdScale*2.*(D.a-C.a);
            
        case 8 % Use no decompression length, free decompression from Woods 1995/Colucci 2014
               % - no dynamic pressure, no velocity limit
            D.pd         = P.pa;
            D            = paramsAtPdFull(D.pd,C,D,P);
            D.u          = C.u + 1./(C.rho_magma.*C.u).*(C.pm-D.pd);
            D.a          = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); % new radius
            D.Ld         = 0;
            
       case 9 % Use no decompression length, free decompression from Woods 1995/Colucci 2014
               % - no dynamic pressure, Decompression velocity limited to Mach 1
            D.pd         = P.pa;
            D            = paramsAtPdFull(D.pd,C,D,P);
            D.u          = D.c; %C.u + 1./(C.rho_magma.*C.u).*(C.pm-D.pd);
            D.a          = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); % new radius
            D.Ld         = 0;
            
    end
    
    if D.Ld<0; D.Ld=0; end % No negative Ld
    
    % Check Ld greater than Zw?
    D.Ld_gt_Zw = D.Ld>=C.Zw;
    
    % Get Ambient pressure at decompression length 
    % If decompression clears water layer, re-equilibrate to
        % atmospheric pressure
    if D.Ld_gt_Zw
        D.pf = P.Patmo;
        
        if D.Ld_gt_Zw && D.pd~=D.pf
            D.pd = P.Patmo;
        end
    else
        D.pf = P.pa - P.g.*P.rho_l.*D.Ld;
    end
    
    % Update params here
    D = paramsAtPdFull(D.pf,C,D,P);
    switch model
        case {2,4,6,7}
            D.u = min([D.c D.u]);
            D.a = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2);
        case {1,3,8}
            D.u = C.u + 1./(C.rho_magma.*C.u).*(C.pm-D.pd);
            D.a = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2); 
        case 5
            D.a = C.a.*(C.pm./D.pd).^(1/2);                      
            D.u = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.a.^2) ); 
    end

end

function D = paramsAtPdFull(P_d,C,D,P)
% Get density, bulk density, gas vol fraction, sound speed, etc etc as f'n of pressure

    equilBubbles = true; % Allow pyroclast bubble gas to equilibrate with ambient? 
                         % (probably not stable in MWI model if false right now)

    W         = waterProps(P_d,C.T);
    
    if equilBubbles   % --> Let bubble gas equilibrate for now, update params
        n_bs0   = sum(P.nsi.*P.ni); % Initial mass frac of bubbles in pyroclasts
        [D.nsi,D.ni,D.Rgsd,D.rhoi,D.pori] = getPSD(P.rho_m, P.phi0, P_d, C.T, 1, P.D, P.phiSz_min, P.phiSz_max);
        n_bs    = sum(D.nsi.*D.ni);
        D.n_0   = ((1-P.n_0)*(1-n_bs0) - 1 + n_bs)./(n_bs - 1);
        D.rho_s   = 1./(sum(D.nsi./D.rhoi));                          % Solids density
        D.Ks      = (D.rho_s.*sum(D.nsi./D.rhoi.*(D.pori./W.Kg + (1-D.pori)./C.Km))).^(-1); % Pyroclast bulk mod
    else
        D.n_0   = C.n_0;
        D.rho_s = C.rho_s;
        D.Ks    = C.Ks;
    end

    D.Kg      = W.Kg;                                               % Adiabatic bulk modulus
    D.rho_g   = W.rho;                                              % Gas density
    D.Kg      = W.Kg;                                               % Gas bulk modulus
    D.rho_d   = 1./( (1-D.n_0)./D.rho_s + D.n_0./D.rho_g);          % Decompression bulk density
    D.chi      = D.rho_d.*D.n_0./D.rho_g;                            % Gas volume fraction
    D.c       = (D.rho_d.*(D.chi./D.Kg + (1-D.chi)./D.Ks)).^(-1/2);   % New decompression sound speed (pseudogas)
%     D.u       = D.c;
%     D.a       = ( C.rho_magma.*C.u.* C.a.^2 ./ (D.rho_d .* D.u) ).^(1/2);
    
end

% function D = paramsAtPd(P_d,C,D)
% % Get density, bulk density, gas vol fraction, sound speed as f'n of pressure
% %     [D.rho_g,D.Kg] = EoS_H2O_2(P_d,C.T);                          % Decompression gas density
%     W            = waterProps(P_d,C.T);
%     D.rho_g      = W.rho;
%     D.Kg         = W.Kg;
%     D.rho_d      = 1./( (1-C.n_0)./C.rho_s + C.n_0./D.rho_g);   % Decompression bulk density
%     D.chi           = D.rho_d.*C.n_0./D.rho_g;                   % Gas volume fraction
%     D.c          = (D.rho_d.*(D.chi./D.Kg + (1-D.chi)./C.Ks)).^(-1/2);  % New decompression sound speed
%     
% end

