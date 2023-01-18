function varargout = Conduit_flow_run


clear
tic


Z0              = 6000;         % Initial depth (m)
dP              = 10e6;         % Initial overpressure (Pa)
T               = 850;          % Temperature (C)
theta           = 145;          % Contact angle between bubbles and crystals
ST_coeff        = 1/4*(2-cos(theta*pi/180))*(1+cos(theta*pi/180))^2;    % Heterogeneous nucleation factor
Q               = 1e8;          % Mass discharge rate (kg/s)
phi_frag        = .75;          % Critical porosity for fragmentation
N0              = 0;            % Initial bubble number density
phi0            = 0;            % Volume fraction of initial exsolved volatiles

conduit_radius  = 37;           % Conduit radius
%==========================================================================
% The run stops when the top boundary condition (pressure = 1 atm or 
% Mach number = 0.99) is reached. 
% One should play around with the conduit radius such that the top boundary
% condition is reached near Z = 0
%==========================================================================




Input  = table(Z0,dP,T,theta,ST_coeff,Q,conduit_radius,phi_frag,N0,phi0);

Output = Conduit_flow_with_nucleation(Input);
    


% save('Results','Output','Input')

disp(['Elapsed time is ' num2str(toc) ' seconds'])

disp(['The run stopped at       Z = ' num2str(Output.Z(end)/1e3) ' km'])

disp(['The final pressure is    P = ' num2str(Output.pm(end)/1e6) ' MPa'])

disp(['The final Mach number is M = ' num2str(Output.M(end))])


if nargout==1
    varargout{1} = Output;
end



end









