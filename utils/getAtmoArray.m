function [A,vars,units] = getAtmoArray(atmtable)
% OUTPUT array: [Press, Alt, Temp Uwind Vwind]

    if isstruct(atmtable)
        vars = fieldnames(atmtable);
        units = atmtable.Units;
        A = zeros(size(atmtable.Pressure,1),5);
        
    elseif istable(atmtable)
        vars = atmtable.Properties.VariableNames;
        atmtable = rmmissing(atmtable);
        units = atmtable.Properties.VariableUnits;
        A = zeros(size(atmtable,1),5);
    end
    
    
    A(:,1) = atmtable.Pressure;
    A(:,2) = atmtable.Altitude;
    A(:,3) = atmtable.Temperature;

    % Just take the vector magnitude of wind
    if ismember('windAbs',vars)
        A(:,4) = atmtable.windAbs;
    else
        wU     = atmtable.MeridionalWindSpeed;
        wV     = atmtable.ZonalWindSpeed;
        A(:,4) = (wU.^2 + wV.^2).^(1/2); 
    end
    A(:,5) = atmtable.relativeHumidity;
    
end