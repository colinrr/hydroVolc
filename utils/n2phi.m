function [phi_i,rho_b] = n2phi(n_i,rho_i)
% Calc mixture mass fractions from volume fractions and densities
%   n_i    Vector or matrix of arbitary n mass fractions
%   rho_i  Vector or matrix of n densities
%
% OUT: phi_i = vector n volume fractions
%      rho_b = bulk density
narginchk(2,2)
assert(all(size(n_i)==size(rho_i)))

rho_b   = sum(n_i(:)./rho_i(:)).^(-1);
phi_i   = rho_b.*n_i./rho_i;

end