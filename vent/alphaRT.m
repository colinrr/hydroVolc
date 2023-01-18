function c_rt = alphaRT(r_0,Tw,p,q,rho_B,u,r)
    % Calculate Rayleigh-Taylor entrainment coefficient for gas jet in
    % water. From Zhang et al 2020.

    % Surface tension params, IAPWS
    Tc = 647.096;
    B  = 0.2358;
    b  = -0.625;
    mu = 1.256;
    tau = 1 - Tw./Tc;
    sigma = B.*tau.^mu .* (1 + b.*tau); % Water surface tension (N/m)
    rho_l = density(p,Tw);

    a = (0.3*u)^2 / (2*pi*r); % Interfacial acceleration
%     a = u^2 / (r);
    c_rt = 4*pi.*r_0.*r./q .* sqrt(2/3 .* rho_B .* (3.*sigma.*rho_l.*a).^(1/2));
end