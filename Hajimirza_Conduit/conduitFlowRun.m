function [cO,cI] = conduitFlowRun(conduitInput)
% Simple wrapper of Conduit flow model to allow sweep-style error handling
% Using Conduit_flow_with_nucleation_V7.
%
% C Rowell, Dec 2023

    if all(isfield(conduitInput, {'atmo','composition','conduit_radius','dP','f0',...
        'mFailTol','N0','pf','pFailTol','phi0','phi_frag','proxy','Q',...
        'rho_melt','ST_coeff','T','theta','vh0','Z0','zFailTol','Zw'}))

        cI = conduitInput;
    else
        cI = getConduitSource(conduitInput);
    end

    try
        cO = Conduit_flow_with_nucleation_V7(cI);
    catch ME
        cO = ConduitOutcome.getErrorOutcomeFields;
        cO.Outcome = ConduitOutcome(ME);
    end

end