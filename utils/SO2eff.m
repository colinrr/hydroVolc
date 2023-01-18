function stratfrac=SO2eff(SO2_height,vent_height,tropo_height)
%SO2_height: the height of SO2 injection in km a.s.l.
%vent_height: the height of the eruptive vent
%tropo_height: the tropopause height at the eruption location in km a.s.l
%stratfrac: the fraction of SO2 injected into the stratosphere


%in the future we could think about adding some function here to account
%for lofting induced by radiation absorption and heating
%SO2_height=cloudloft(SO2_height,otherinputs);


%SO2 injection profile parameterized using 3D model outputs, see Figure S2
%in Aubry et al. (GRL 2019)
SO2prof = @(z) exp(-(z-SO2_height).^2./(0.108*(SO2_height-vent_height)).^2);

%calculate the fraction above the tropopause
stratfrac=integral(SO2prof,tropo_height,Inf)./integral(SO2prof,-Inf,Inf);


end