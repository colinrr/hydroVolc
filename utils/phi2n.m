function [n_i,rho_b] = phi2n(phi_i,rho_i)
% Calc mixture mass fractions from volume fractions and densities
%   phi_i  Vector or matrix of arbitary n volume fractions
%   rho_i  Matching vector or matrix of n densities
%
% OUT: n_i = vector n mass fractions
%      rho_b = bulk density
narginchk(2,2)
assert(all(size(phi_i)==size(rho_i)))

rho_b   = dot(phi_i,rho_i);
n_i     = rho_i.*phi_i./rho_b;

end