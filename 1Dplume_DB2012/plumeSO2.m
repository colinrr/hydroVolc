function [F, m_SO2] = plumeSO2(pI,pO)
% Function to calculate in-plume SO2 scavenging and stratospheric delivery
% efficiency.
%  --> Need SO2_0 in pI

    Mso2 = 0.064066;    % SO2 molar mass, kg/mol
    Rso2 = 8.314/Mso2;  % SO2 gas constant
    gamma = 6.2e-8;     % Uptake coefficient

    % Particle Specific Surface Area - geometric
    SSAi  = 3./(pI.Rgsd.*pI.rhoi)'.*pO.mi; % by grain size [m^2/kg]
    SSA   = sum(SSAi,2); % Specific surface area of bulk particles [m^2/kg]
    SA_B = SSA.*pO.m_s./pO.m.*pO.rho_B; % Bulk mixture volumetric SA [m^2/m^3]
    
    cso2 = sqrt(3*Rso2*pO.theta); % Molecular velocity... 

    % SSA BET
%     SSA  = 500; % m^2/kg
%     gamma = 1e-8;

    
%     SSAi = 3./(pI.Rgsd.*pI.rhoi)'.*pO.mi;
%     rho = 2400;
%     dp   = 10e-6; % avg particle diameter
%     SAp  = 4*pi*(dp/2)^2;
    % Vp   = 4/3*pi*
%     SSA = 3/(dp/2)/rho;
%     gamma = 1e-3;

%     rhop = .002; % kg/m3 particle loading

    % Something about this equation cannot be right. Must be [SO2] in here
    % somewhere?
    tau = 4./(gamma.*cso2.*SA_B);
    
    % dSO2/dz = -SO2
    % SA = SSA*pO.m_s;
    % dSA/dz = -SED - dSO2*(SA/SO2)

    % SO2 delivery efficiency from gaussian height distribution
    F = SO2eff(pO.hb+pI.vh0, pI.vh0, pO.atmo.ztropo);
%     F = F.*SO2(end)./pI.SO2_0;
end

function stratfrac=SO2eff(SO2_height,vent_height,tropo_height)
    %SO2_height:    the height of SO2 injection in km a.s.l.
    %vent_height:   the height of the eruptive vent
    %tropo_height:  the tropopause height at the eruption location in km a.s.l
    %stratfrac:     the fraction of SO2 injected into the stratosphere

    %in the future we could think about adding some function here to account
    %for lofting induced by radiation absorption and heating
    %SO2_height=cloudloft(SO2_height,otherinputs);

    %SO2 injection profile parameterized using 3D model outputs, see Figure S2
    %in Aubry et al. (GRL 2019)
    SO2prof = @(z) exp(-(z-SO2_height).^2/(0.108*(SO2_height-vent_height))^2);

    %calculate the fraction above the tropopause
    stratfrac=integral(SO2prof,tropo_height,Inf)/integral(SO2prof,-Inf,Inf);


end