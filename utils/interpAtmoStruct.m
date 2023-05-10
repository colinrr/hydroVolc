function A = interpAtmoStruct(atmo,z)
% input atmo struct, fields assumed in SI units

A = atmo;
Z = atmo.Altitude;
interpMethod = 'pchip';

if isfield(atmo,'Altitude')
    A.Altitude = z;
end
if isfield(atmo,'Pressure')
    A.Pressure = 10.^interp1(Z,log10(atmo.Pressure),real(z),interpMethod,'extrap');
end
if isfield(atmo,'Temperature')
    A.Temperature = interp1([Z;6e4;8.5e4],[atmo.Temperature;273;173],real(z),interpMethod,'extrap');
end
if isfield(atmo,'MeridionalWindSpeed')
    A.MeridionalWindSpeed = twoSidedExtrap(Z,atmo.MeridionalWindSpeed,real(z),interpMethod,atmo.MeridionalWindSpeed(end));
%     A.MeridionalWindSpeed = interp1(Z,atmo.MeridionalWindSpeed,real(z),interpMethod,atmo.MeridionalWindSpeed(end));
end
if isfield(atmo,'ZonalWindSpeed')
    A.ZonalWindSpeed = twoSidedExtrap(Z,atmo.ZonalWindSpeed,real(z),interpMethod,atmo.ZonalWindSpeed(end));
%     A.ZonalWindSpeed = interp1(Z,atmo.ZonalWindSpeed,real(z),interpMethod,atmo.ZonalWindSpeed(end));
end
if isfield(atmo,'windAbs')
    A.windAbs = twoSidedExtrap(Z,atmo.windAbs,real(z),interpMethod,atmo.windAbs(end));
%     A.windAbs = interp1(Z,atmo.windAbs,real(z),interpMethod,atmo.windAbs(end));
end
if isfield(atmo,'relativeHumidity')
    A.relativeHumidity = interp1(Z,atmo.relativeHumidity,real(z),interpMethod,'extrap');
end

end

function y = twoSidedExtrap(Z,Y,z,interpMethod,extVal)
% Extrapolate with specific value only on the high side
% Z,Y: input vectors
% z: output z vector
% extVal = extrap value for interp function
y = zeros(size(z));
zCheck = z<min(Z);
    if any(zCheck)
        y(zCheck) = interp1(Z,Y,real(z(zCheck)),interpMethod,'extrap');
    end
    y(~zCheck) = interp1(Z,Y,real(z(~zCheck)),interpMethod,extVal);
end