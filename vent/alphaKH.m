function c_kh = alphaKH(r_0,q_0,rho_B,u,r)
    % Calculate Kelvin-Helmholtz entrainment coefficient for gas jet in
    % water. From Zhang et al 2020.
    
    c_kh = 2*pi/sqrt(3) .* (2.*r_0) ./ q_0 .* rho_B .* u .* r;

end