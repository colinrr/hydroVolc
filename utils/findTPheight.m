function htropo=findTPheight(alt,temp)
%alt=altitude in km a.s.l.
%temp=temperature in degree Celsius or Kelvin
%both alt and temp should be 1D arrays of the same size
%htropo will be the thermal tropopause height (WMO definition) in km a.s.l.



%interpolate temperature profile with 50m resolution between 5 and 25km
%(assumes tropopause height within these bounds)
zref=5:0.05:25;
temp=interp1(alt,temp,zref,'pchip');

%calculate temp gradient
dtdz=diff(temp)./diff(zref);
%shift altitude to mid-point for consistency with gradient calculation
zref=zref(1:end-1)+0.025;

%find where temperature gradient exceeds 2degree/km
indgrad=find(dtdz>=-2);

%==========================================================================
%the whole while loop thing below finds the first height at which it exceeds 2deg/km
%for at least 2 km (40*50m)
itest=1;
ztr=0;
while ztr==0 & itest<=length(indgrad)-39
    dtdztest=dtdz(indgrad(itest):indgrad(itest)+39);
    if isempty(dtdztest(dtdztest<-2))
        ztr=1;
    else
        itest=itest+1;
    end
end
htropo=zref(indgrad(itest));

end