function A = interpAtmoArray(atmo,z)
% input atmo struct, fields assumed in SI units

A = zeros(length(z),5);
Z = atmo(:,2);
P = atmo(:,1);
T = atmo(:,3);
V = atmo(:,4);
RH = atmo(:,5);
interpMethod = 'pchip';

    A(:,2) = z;
    A(:,1) = 10.^interp1(Z,log10(P),real(z),interpMethod,'extrap');

    A(:,3) = interp1([Z;6e4;8.5e4],[T;273;173],real(z),interpMethod,'extrap');

    A(:,4) = twoSidedExtrap(Z,V,real(z),interpMethod,V(end));

    A(:,5) = interp1(Z,RH,real(z),interpMethod,'extrap');


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