function phi_s0 = sourceParticleVolFrac(dat)
% Calculate particle volume fraction at plume source.

N = numel(dat);
phi_s0 = nan(size(dat));
for ni = 1:N
    if ~dat(ni).wO.failedPlume
        % Mass fractions for all particle phases at source relative to total
        n_si    = (dat(ni).pO.m_si(1,:)./dat(5).pO.m(1))'; % Total particle mass fractions
        rho_si  = dat(ni).pI.rhoi;  % Various particle densities

        % Add gas and liquid props
        n_all   = [dat(ni).pI.n_0; n_si];
        rho_all = [dat(ni).pI.rho_g0; rho_si];

        phi_all = n2phi(n_all,rho_all); % All volume fractions
        phi_s0(ni)  = sum(phi_all(2:end)); % Total particle volume fraction only
    end
end


end
