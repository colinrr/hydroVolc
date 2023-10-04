function [Zpass,UnderPressurePass,ChokePass,PressBalancePass,FragCheck,FlareCheck,valid,report] = checkConduitResult(cO,Zthresh,Mthresh,Pthresh)
% [Z0pass,UnderPressurePass,ChokePass,PressBalancePass,FragCheck,FlareCheck,valid] = checkConduitResult(cO,Zthresh,Mthresh,Pthresh)
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
%   Z0pass              1) TRUE if upper-most z value is less than Zthresh
%                       ie Conduit run reaches surface. Typical Zthresh is 
%                       within a conduit diameter or so
%
%   UnderPressurePass  2) TRUE if (P_v + rho*u^2/2)/P > (1- Pthresh) : Conduit is not
%                       underpressured (ie likely to collapse).
%
%  AND EITHER OF
% ChokePass             3) TRUE if Mach # < Mthresh with P_v/P_a > (1-Pthresh) :
%                        Conduit is choked and either pressure balanced or overpressured
% PressBalancePass    4) TRUE if Mach # < Mthresh and |dP/P| < (1 +/- Pthresh)
%                       Conduit is not choked but is pressure balanced
%
% Additional diagnostics:
% FragCheck    5) Did fragmentation occur?
% FlareCheck   6) Does conduit flare? ie vent radius > initial radius
%
% valid = true/false for successful run based on the above conditions

if nargin<5 || isempty(verbose)
    verbose = false;
end

% Did fragmentation occur?
FragCheck   = logical(cO.Par.frag);

Zpass   = cO.Z(end) < Zthresh;                                                     % Reaches surface
UnderPressurePass  = ( (cO.pm(end)+cO.U(end).^2*cO.rho_magma(end)/2)/cO.Par.pf ) > (1-Pthresh); % Not underpressured
PressBalancePass   = and( (cO.pm(end)/cO.Par.pf) < (1+Pthresh), UnderPressurePass );              % Pressure balanced
ChokePass   = and( cO.M(end)>Mthresh,cO.M(end)<1.1 );                           % Choked - upper limit a bit arbitrary here as it isn't really used

% Fcheck  = cO.a(end) > cO.a(1);
FlareCheck  = abs((cO.a(end) - cO.a(1))./cO.a(1)) > .005;

% Non - fragmented case corresponds to and(Zpass,UnderPressurePass) && PressBalancePass
%   --> unless we allow flow to not reach surface...

% Fragmented case
valid   = and(Zpass,UnderPressurePass) && or(PressBalancePass,ChokePass);

report = sprintf('R= %.5f, Z ~ 0? %i, Frag? %i, UnderP.? %i, P bal.? %i, Mach#? %i, Flare? %i, Valid? %i\n',...
        cO.Par.a,Zpass,FragCheck,UnderPressurePass,PressBalancePass,ChokePass,FlareCheck,valid);
 

end