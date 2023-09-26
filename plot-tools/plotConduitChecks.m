
% Testing some params to check for and plot valid conduit solutions

dd = D(1);

% Checks: Frag, choke, flare
% Minimization plots needed:
    % Z0
    % P - split into OP/UP?
    
% x axes options: varied param, or run order   

for di = 1:length(D)
    dd = D(di);
    
    % Get thresholds for pressure, mach number, z
    Z_thresh_m(di) = dd.cI.ZfailScale*dd.cI.conduit_radius;
%     UP_thresh_pa   = dd.
    
    
    [Z0pass,UnderPressurePass,ChokePass,PressBalancePass,FragCheck,FlareCheck,valid,report] = ...
        checkConduitResult(dd.cO,dd.cI.ZfailScale*dd.cI.conduit_radius,dd.cI.Mfailthresh,dd.cI.Pfailthresh);

end



figure
subplot()

