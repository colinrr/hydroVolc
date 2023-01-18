function dat = reRunMCsingle(cI,pI,MC,randPars,idx)
% Given a set of input params for a monte carlo simulation, re-run a single
% simulation using the parameter set corresponding to index "idx"

cff = fieldnames(MC.cI);
pff = fieldnames(MC.pI);

for ci = 1:length(cff)
    
    if strcmp(cff{ci},'logQ')
%         cI.Q = randPars(idx,'Q');
        cff{ci} = 'Q';
    elseif strcmp(cff{ci},'a_var')
%         cI.conduit_radius = randPars(idx,'conduit_radius');
        cff{ci} = 'conduit_radius';
    end
    cI.(cff{ci}) = randPars{idx,cff{ci}};
%     end
end

for pi = 1:length(pff)
    pI.(pff{pi}) = randPars{idx,pff{pi}};
end


dat = runCoupledModel(cI,pI);

end