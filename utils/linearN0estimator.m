function N0 = linearN0estimator(Q,P)
% This is a very simple multiple linear regression fit for initial BND
% after a single nucleation event, based on output of the conduit model. It
% is intended as a simple estimator for an appropriate value of N0 given a
% non-zero initial exsolved fraction (phi0). 
% The two fit variables are Q (mass discharge rate) and pf...
%  - See bnd_h2o_MLRFits.m and N_vs_phi_checks.m for origin.
%  - A more complete MLR using Q, P0 (dP + pf + rho*g*Z0), composition could
%    be useful in future.

p00 = 10.076302756178352;
pQ  = 0.514654850750388;
pP  = -5.814103285783626e-08;

lQ = log10(Q);

lN = p00 + pQ.*lQ + pP.*P;

N0 = 10.^lN;

end