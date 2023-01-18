function [pI,dO,wO] = con2plume(cI, cO, pI)
% plumeSource = con2plume(conduitI, conduitO);
% Convert conduit output into input for plume via...
% 1 Vent Decompression
% 2 External Water Interaction
%
%   cI = table from getConduitSource
%   cO = output struct from Conduit_flow_with_nucleation
%   pI = manual input params for plume model (struct, see getPlumeSource).
%        MOST will be overwritten here from conduit output.
%           D, phiSz_min, phiSz_max

    narginchk(2,6)
    if nargin<3
        pI = [];
    end
    
    if cO.Par.frag

        % No GSD from conduit, for now:
        %--- Retrieve params from Sahand's conduit model at vent ---
%         [rho_in,Kg_in] = EoS_H2O_2(cO.pg(end),cO.Par.T);
        Win = waterProps(cO.pg(end),cO.Par.T);
        rho_in = Win.rho;
        Kg_in  = Win.Kg;
        vg_in = 4/3 * pi * cO.M3(end);
        porosity  = vg_in./(1+vg_in);  % true BUBBLE porosity (fraction of pyroclast volume)
        
        % Bulk volume fractions
        phi_s  = 1 - cO.porosity(end);  % Solids (no bubbles) total volume fraction
        phi_b  = (phi_s*porosity)/(1-porosity);       % Bubbles total volume fraction
        phi_f  = 1 - phi_b - phi_s;     % Free gas total volume fraction

        m_s    = cO.Par.rho_melt * phi_s / cO.rho_magma(end) * cO.Par.Q;    % Mass of solids
        m_vf   = cO.rho_g(end) *  phi_f / cO.rho_magma(end) * cO.Par.Q;     % Free gas mass
        m_vb = rho_in * phi_b / cO.rho_magma(end) * cO.Par.Q;               % Bubble gas mass
        m_v  = m_vb + m_vf;                                                 % Total exsolved gas mass
        n_vf = m_vf/cO.Par.Q;                                               % Free gas mass fraction
        ng = m_v/cO.Par.Q;                                                  % Total exsolved gas mass fraction
        
        % Initiate plume source with inherited conduit params
        pI.T0        = cO.Par.T;
        pI.rho_m     = cO.Par.rho_melt;
        pI.atmo      = cI.atmo;
        pI.vh0       = cI.vh0;
        if ~isfield(pI,'T_g')
            pI.T_g       = 785.5 - 83.48*log(cO.Cm(end)*100); % Dingwell 1998 rhyolite glass transition, fit corrected
        end
        if ~isempty(pI)
            pI = getPlumeSource(pI);
        else
            pI = getPlumeSource;
        end
        
        % Generate grain size distribution
        [pI.nsi,pI.ni,pI.Rgsd,pI.rhoi,pI.pori] = getPSD(cO.Par.rho_melt, porosity, cO.pg(end), cO.Par.T, 1, pI.D, pI.phiSz_min, pI.phiSz_max);
        
        % -------- NEW params following 'application' of PSD ----------
        pI.n_0 = ng - (1-ng).*sum(pI.ni.*pI.nsi);                            % New free gas mass fraction
        rho_bs2 = 1./(sum(pI.nsi./pI.rhoi));                                 % New bulk solids
        rho_g2 = pI.n_0.*cO.rho_g(end).*phi_f./n_vf.*...                     % New free gas density
            (1 - cO.rho_g(end).*phi_f.*(1-pI.n_0)./(n_vf.*rho_bs2)).^-1;
%         [Pv,~] = EoS_H2O(rho_g2,cO.Par.T);                                 % New vent gas pressure
        Pv     = pressure(rho_g2,cO.Par.T);
        rho_b2 = 1./(pI.n_0/rho_g2 + (1-pI.n_0)./rho_bs2);                   % New bulk density
        
        % NOTE: after this stage of PSD initialization, calculated "solids"
        % mass should go UP a bit, because it now includes mass of both
        % bubbles and melt. Previously it was melt only, and all exsolved
        % gass.
        
        %% -------- Run vent decompression ----------
        % (1) Estimate Ld scale and select decompression regime
        % (2) Get decompression params: P, u, rho, etc
        % (3a) If Ld < Zw: Run MWI to Zw
        % (4) Go to plume model
        

        % Decompression input
        dI.Zw        = cI.Zw;
        dI.rho_magma = cO.rho_magma(end);
        dI.rho_g     = rho_g2;
        dI.rho_s     = rho_bs2;
%         dI.C         = cO.C(end); % Approximate
        dI.u         = cO.U(end);
        dI.c         = cO.C(end);
        dI.n_0       = pI.n_0;
        dI.a         = cO.a(end);
        dI.pm        = Pv;
        dI.pf        = cO.Par.pf;
        dI.T         = cO.Par.T;
        dI.Km        = cO.Par.K_melt;
        dI.Ks        = (rho_bs2.*sum(pI.nsi./pI.rhoi.*(pI.pori./Kg_in + (1-pI.pori)./cO.Par.K_melt))).^(-1);
        
        [dO,pI] = getLd(dI,pI); % Vent decompression
        pI.vh0       = cI.vh0 + dO.Ld; %cI.Zw;
    

%% Magma-water interaction model

        if cI.Zw>pI.Ld
            % Choose MWI model
            switch pI.mwiModel
                case 1
                    [pI,wO] = MWIv1(cI,cO,pI,dO);
                case 2
                    [pI,wO] = MWIv2(cI,cO,pI);
                case 3
                    [pI,wO] = MWIv3(cI,cO,pI,dO); % IN DEV...
                otherwise
                    error('mwiModel version not recognized.')
            end
        else
            wO.failedPlume = false;
        end

        % Run plume
%         pI = getPlumeSource(pI);

        %% Mass continuity check

        % Conduit
        % phi = conduitO.porosity(end);
        % m_n_conduit = pi.*dat.conduitO.a(end).^2.*dat.conduitO.U(end).*(phi*dat.conduitO.rho_g(end) + (1-phi)*dat.conduitO.Par.rho_melt);
%         m_n_conduit = pi.*cO.a(end).^2.*cO.U(end).*(cO.rho_magma(end));

        % Plume
%         m_0_plume = pi.*pI.r_0.^2.*pI.u_0.*rhoB_out;

        
%         assert(round(m_n_conduit)==round(m_0_plume),sprintf('Conduit exit mass flux:\t%.2f kg/s\nPlume source mass flux:\t%.2f kg/s\n',m_n_conduit,m_0_plume))

        % fprintf('',m_0_plume)

    else
        pI.T0       = NaN;
        pI.n_0       = NaN; 
        pI.u_0       = NaN;
        pI.r_0       = NaN;
        pI.rho_s    = NaN;
        pI.atmo     = NaN;
        pI.vh0      = NaN;

        % For later comparison
        pI.rho_B0   = NaN;
        pI.rho_g0   = NaN;
        

    end


end
