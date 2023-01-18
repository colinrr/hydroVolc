function varargout = renameOutputVars(qSet,qA,plmIn)


    if ischar(qSet)
        iDat = qSet;
        load(qSet,'qSet');
    end
    if ischar(qA)
        iSumm = qA;
        load(qA,'qA')
    end
    assert(and(isstruct(qSet),isstruct(qA)),'Model output files did not load properly')

%%

    qSetrm.cI = {'atmo'};

    % Fields to rename
    plmInN = {'Tw0','Twe'};
    qSetN.cI = {'pf','pe'
                'phi0','chiI'
                'phi_frag','chi_frag'
                'Zw','Ze'};

    qSetN.Par = {'pf','pe'
                 'mu','eta'};

    qSetN.cO = {'D','Dw'
                'pg','pv'
                'mu','eta'
                'rho_g','rho_v'};


    qSetN.dO = {'pd','pd0'
                'ni','nbi'
                'Kg','Kv'
                'rho_g','rho_v'};

    qSetN.wO = {'phi_pdf','nsif'
                'm_s','q_s'
                'psi','qu'
                'm','q'
                'r_mix','a_mix'
                'm_v','q_v'
                'm_l','q_l'
                'r','a'};

    qSetN.pI = {'beta','b'
                'ni','nbi'
                'phi0','chi0'
                'r_0','a_0'
                'Tw0','Twe'};

    qSetN.pO = {'hm','zm'
                'rm','am'
                'hb','zb'
                'm_0','q_0'
                'm','q'
                'm_d','q_d'
                'm_l','q_l'
                'm_v','q_v'
                'm_si','q_si'
                'm_s','q_s'
                'psi','qu'
                'Q','E'
                'r','a'
                'theta','T'};

    qAN.cI = {'phi_frag','chi_frag'
        'Zw','Ze'
        'pf','pe'};

    qAN.cO = {'pf','pe'
              'mu_fr','eta_fr'
              'pg_fr','pv_fr'};

    qAN.dO = {'pd','pd0'};

    qAN.wO = {'m','q'
              'm_s','q_s'
              'm_l','q_l'
              'm_v','q_v'};

    qAN.pI = {'r_0','a_0'};

    qAN.pO = {'hm','zm'
              'hb','zb'
              'rm','am'
              'm','q'
              'm_d','q_d'
              'm_l','q_l'
              'm_v','q_v'
              'theta','T'
              'm_0','q_0'
              'SSA_hm','SSA_zm'
              'SSA_hb','SSA_zb'
              'SAv_hb','SAv_zb'
              'm_w_tp','q_w_tp'
              'm_s_tp','q_s_tp'};

%% Rename qSet vars
    pl = size(plmInN,1);
    for fi = 1:pl
        plmIn = renameStructField(plmIn,plmInN{fi,1},plmInN{fi,2});
    end
    
    
    fn = fieldnames(qSetN);
    for ff = 1:length(fn)
        fi = fn{ff};

        varSet = qSetN.(fi);
        nv = size(varSet,1);

        for vi = 1:nv
            for nn = 1:numel(qSet)
                switch fi
                    case 'Par'
                        if ~isempty(qSet(nn).cO)
                            try
                                qSet(nn).cO.(fi) = renameStructField(qSet(nn).cO.(fi),varSet{vi,1},varSet{vi,2});
                            catch
                                fprintf('%s %s %s %i',fi,varSet{vi,1},varSet{vi,2},nn)
                                flargh
                            end
                        end
                    otherwise
                        if ~isempty(qSet(nn).(fi))
                            try
                                qSet(nn).(fi) = renameStructField(qSet(nn).(fi),varSet{vi,1},varSet{vi,2});
                            catch
                                fprintf('%s %s %s %i',fi,varSet{vi,1},varSet{vi,2},nn)
                                flargh
                            end
                        end
                end
            end
        end    
    end
   
    %% Rename qA vars
    
    fn = fieldnames(qAN);
    for ff = 1:length(fn)
        fi = fn{ff};
        
        varSet = qAN.(fi);
        nv = size(varSet,1);
        
        for vi = 1:nv
            qA.(fi) = renameStructField(qA.(fi),varSet{vi,1},varSet{vi,2});
        end    
    end
    qA = renameStructField(qA,'Zw','Ze');

%% delete vars
rf = fieldnames(qSetrm);
for ri=1:length(rf)
    rv = qSetrm.(rf{ri});
    for rvi = 1:length(rv)
        for nn = 1:numel(qSet)
            if ~isempty(qSet(nn).(rf{ri}))
                qSet(nn).(rf{ri}) = rmfield(qSet(nn).(rf{ri}),rv{rvi});
            end
        end
    end
end

%%
    if nargout >2
        varargout{3} = plmIn;
    end
    if nargout >1
        varargout{2} = qA;
    end
    if nargout >0
        varargout{1} = qSet;
    end
end