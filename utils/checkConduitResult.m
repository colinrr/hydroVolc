function [Zpass,UPpass,Mpass,Ppass,Fcheck,valid] = checkConduitResult(cO,Zthresh,Mthresh,Pthresh)
% [Zpass,UPpass,Mpass,Ppass,Fcheck,valid] = checkConduitResult(cO,Zthresh,Mthresh,Pthresh)
% Check conduit output for validity
% IN:   cO      = conduit model output struct
%       Zthresh = maximum conduit "top depth" (m). 
%                   Default = 2*r0 (1 conduit diameter)
%       Mthresh = Check Mach number great than Mthresh.
%                   Default = 0.95;
%       Pthresh = Check abs((P_vent - P_surface)/P_surface)< Pthresh. 
%                 ie check if vent and surface pressures differ by less
%                 than this fraction.
%                   Default = 0.1;
%
% PASS CONDITIONS ((ALWAYS REQ'D):
%   Zpass  1) TRUE if upper-most z value is less than Zthresh
%                       ie Conduit run reaches surface. Typical Zthresh is 
%                       within a conduit diameter or so
%
%   UPpass  2) TRUE if (P_v + rho*u^2/2)/P > (1- Pthresh) : Conduit is not
%               underpressured (ie likely to collapse).
%
%  AND EITHER OF
% Mpass    3) TRUE if Mach # < Mthresh with P_v/P_a > (1-Pthresh) :
%               Conduit is choked and either pressure balanced or overpressured
% Ppass    4) TRUE if Mach # < Mthresh, |dP/P| < (1-Pthresh)
%               Conduit is not choked but is pressure balanced
%
% Additional diagnostic:
% Fcheck   5) Does conduit flare? ie vent radius > initial radius
%
% valid = true/false for successful run based on the above conditions

Zpass   = cO.Z(end) < Zthresh;                                              % Reaches surface
UPpass  = ( (cO.pm(end)+cO.U(end).^2*cO.rho_magma(end)/2)/cO.Par.pf ) > (1-Pthresh); % Not underpressured
Ppass   = and( (cO.pm(end)/cO.Par.pf) < (1+Pthresh), UPpass );              % Pressure balanced
Mpass   = and( cO.M(end)>Mthresh,cO.M(end)<1.1 );                           % Choked

% Fcheck  = cO.a(end) > cO.a(1);
Fcheck  = abs((cO.a(end) - cO.a(1))./cO.a(1)) > .005;

valid   = and(Zpass,UPpass) && or(Ppass,Mpass);

end