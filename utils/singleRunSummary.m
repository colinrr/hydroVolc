function T = singleRunSummary(dat)
% Prints a summary table for a single model run
% d = coupledModel_singleRun output struct
%
% CRowell Mar 2021

c = 'cExit';
p = 'pSource';
T = table;

T{'MER',{c p}} =    [pi.*dat.conduitO.a(end).^2.*dat.conduitO.U(end).*(dat.conduitO.rho_magma(end))
                    dat.plumeO.m_0]';

T{'U',{c p}}    =   [dat.conduitO.U(end) 
                    dat.plumeI.u_0]';

T{'r',{c p}}    =   [dat.conduitO.a(end)
                    dat.plumeI.r_0]';

T{'rhoB',{c p}} =   [dat.conduitO.rho_magma(end)
                    dat.plumeI.rho_B0]';

T{'Pm',{c p}}    =  [dat.conduitO.pm(end)
                    dat.plumeO.P_0]';

T{'Zfrag',{c p}} =   [dat.conduitO.Par.Zf
                      NaN]';

T{'Ztop',{c p}}  =   [dat.conduitO.Z(end)
                      dat.plumeO.hm]';

T{'NBL',{c p}}   =  [NaN
                    dat.plumeO.hb]';

T{'Clps',{c p}}  =  [NaN
                    dat.plumeO.collapse]';


% T{} = ;
end