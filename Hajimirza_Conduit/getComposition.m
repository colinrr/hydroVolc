function composition = getComposition(varargin)
% composition = getComposition(varargin)
    % Magma composition model used to compute viscosity. Currently built
    % using the emprical fit model of Hui and Zhang, 2007.
    % A switch to use the model of Giordano, Russell, and Dingwell 2008 is
    % an obvious dev step.
    %
    % Compositional parameters should be entered as WEIGHT FRACTIONS.
    %
    % Handy list of input parameters and default values:
% composition = getComposition(...
%     'SiO2', 76.53e-2, ...
%     'TiO2', .06e-2, ...
%     'Al2O3', 13.01e-2, ...
%     'FeO',  .79e-2, ...
%     'Fe2O3', 0, ...
%     'MnO',  .08e-2, ...
%     'MgO',  .02e-2, ...
%     'CaO',  .74e-2, ...
%     'Na2O', 3.87e-2, ...
%     'K2O',  4.91e-2 ...
%     );

%% DETAULT INPUT PARAMETERS

SiO2    = 76.53e-2;
TiO2    = .06e-2;
Al2O3   = 13.01e-2;
FeO     = .79e-2;
Fe2O3   = 0;  % The Hui and Zhang model considers all Fe as FeO, and any input Fe2O3 will be converted to total FeO in the viscosity model
MnO     = .08e-2;
MgO     = .02e-2;
CaO     = .74e-2;
Na2O    = 3.87e-2;
K2O     = 4.91e-2;

%% Parse
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    
    addParameter(p,'SiO2',SiO2)
    addParameter(p,'TiO2',TiO2)
    addParameter(p,'Al2O3',Al2O3)
    addParameter(p,'FeO',FeO)
    addParameter(p,'Fe2O3',Fe2O3)
    addParameter(p,'MnO',MnO)
    addParameter(p,'MgO',MgO)
    addParameter(p,'CaO',CaO)
    addParameter(p,'Na2O',Na2O)
    addParameter(p,'K2O',K2O)

    parse(p,varargin{:})
    composition = p.Results;

    

end