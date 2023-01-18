% Convert atmo csv to ready .mat format

iDir = '/Users/crrowell/code/research-projects/glaciovolc/glaciovolc-dev/1Dplume_DB2012';


% iFile = 'atm_ERAreanalysis_Grimsvotn2011_01.csv';
% iFile = 'atm_ERAreanalysis_Shinmoedake2011_02.csv';
iFile = 'atm_ERAreanalysis_Tungarahua2014_01.csv';


%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Altitude", "Temperature", "Pressure", "Zonalwindspeed", "Meridionalwindspeed", "Relativehumidity"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];
VariableUnits = ["m a.s.l.", "K", "Pa", "m/s", "m/s", "%"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
atmprofile = readtable(fullfile(iDir,iFile), opts);
atmprofile.Properties.VariableUnits = VariableUnits;
% Clear temporary variables
clear opts

%% 
[~,oName,~] = fileparts(iFile);
save(fullfile(iDir,oName),'atmprofile','VariableUnits')