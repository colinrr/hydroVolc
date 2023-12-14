classdef (ConstructOnLoad = true) ConduitOutcome
%   ConduitOutcome
%   ConduitOutcome(code)  : Code = integer
%   ConduitOutcome(name)  : Name = string
%   ConduitOutcome(cI,cO) : cI, cO are I/O structs from Conduit model
%   ConduitOutcome(ME)    : ME is MatlabException from failed model run
%
%   Produces and object containing code ID and explanatory strings for 
%   outcome of a conduit model run. 
%   See method ConduitOutcome.printTable for a list of names, codes, and
%   descriptions.

    properties
        Name (1,1) string                   = "nullResult"; %
        Code (1,1) double {mustBeInteger}   = 0; % 
        Label (1,1) string                  = "Code 0 - no result";
        Msg  char                           = 'Code 0 - default constructor, no result'; %
        Exception  % Exception funtionality not yet in
        ZFailTol
        MFailTol
        PFailTol
        Failed                      = false;
        DepthFlag logical           = false;
        NotUnderPressured logical   = false;
        PressureBalanced logical    = false;
        Valid  logical              = false;
        Frag  logical               = false;
        PhiCritReached  logical     = false;
        Choked  logical             = false;
        Flared  logical             = false;
    end
    
    methods
        function obj = ConduitOutcome(varargin)
                
            
            getFromCode = false;
            getFromName = true;
%             getTols     = true;
%             Tols = [];
            

            if nargin>=3
                assert(islogical(varargin{3}),'Third argument is verbose: true/false.')
                verbose = varargin{3};
            else
                verbose = false;
            end

            if nargin == 0
                obj = ConduitOutcome.getDefaultTolerances(obj);
            end
            
            if nargin==1
                if isa(varargin{1},'MException')
                    obj = errorOutcome(obj,varargin{1});
                
                elseif isnumeric(varargin{1}) && isscalar(varargin{1})
                    obj.Code = varargin{1};
%                     getFromCode = true;
                    referenceTable = ConduitOutcome.getTable(obj);
                    obj = getRowByCode(obj,referenceTable);
                    
                elseif ischar(varargin{1}) || isstring(varargin{1})
                    obj.Name = varargin{1};
                    referenceTable = ConduitOutcome.getTable(obj);
                    obj = getRowByName(obj,referenceTable);
                    
                end
            
            elseif nargin >= 2
                cI = varargin{1}; cO = varargin{2}; % These can be classes eventually
                assert(all([isstruct(cI) isstruct(cO) isfield(cI,'conduit_radius') isfield(cO,'rho_magma')]),...
                    'Two arguments must be input and output structs of conduit model.')
            
                referenceTable = ConduitOutcome.getTable(obj);

                [obj,~] = ConduitOutcome.checkConduitTolerances(cI,cO,obj);

                obj.Name = ConduitOutcome.interpretConduitResult(obj,verbose);

                if getFromCode
                    obj = getRowByCode(obj,referenceTable);

                elseif getFromName
                    obj = getRowByName(obj,referenceTable);

                end
            end
            
        end
        
        function obj = errorOutcome(obj,ME)
            knownIdentifiers = {'physicalBoundsError','IntegrationTolNotMet'};
            referenceTable = ConduitOutcome.getTable(obj);

            if contains(ME.identifier, knownIdentifiers)
                switch ME.identifier
                    case 'MATLAB:ode15s:IntegrationTolNotMet'
                        obj.Name   = 'failedIntegration';
                    case {'findPb:physicalBoundsError','after_fragmentation:physicalBoundsError'}
                        obj.Name   = 'failedPhysics';
                end
                
            else
                obj.Name   = 'failedUnmappedError';
            end

            obj = getRowByName(obj,referenceTable);
            obj.Label  = ME.identifier;
            obj.Msg    = ME.message;
            obj.Exception = ME;
            obj = obj.getDefaultTolerances(obj);
            obj.Failed = true;
        end
        
        function obj = getRowByCode(obj, T)
            assert(ismember(obj.Code,T.Code),'Conduit: Outcome Code not recognized.')
            [~,idx] = ismember(obj.Code,T.Code);
            
            obj.Name = T.Properties.RowNames{idx};
            obj.Label = T{idx,'Label'};
            obj.Msg  = T{idx, 'Message'};
        end
        
        function obj = getRowByName(obj, T)
            assert(ismember(obj.Name,T.Properties.RowNames),'Conduit: Outcome Name not recognized.')
%             [~,idx] = ismember(obj.Name,T.Properties.RowNames);
            
            obj.Code = T{obj.Name,'Code'};
            obj.Label = T{obj.Name,'Label'};
            obj.Msg  = T{obj.Name, 'Message'};
        end
        
        function repString = reportString(obj) %, cI)
%             if isempty(cI)
                repString = sprintf('%s: Z ~ 0: %i, Frag: %i, UnderP.: %i, P bal.: %i, Mach#: %i, Flare: %i, Valid: %i',...
                    obj.Label,obj.DepthFlag,obj.Frag,obj.NotUnderPressured,obj.PressureBalanced,obj.Choked,obj.Flared);
%             else
%                 repString = sprintf('%%s: Z ~ 0: %i, Frag: %i, UnderP.: %i, P bal.: %i, Mach#: %i, Flare: %i, Valid: %i\n',...
%                     obj.Label,obj.DepthFlag,obj.Frag,obj.NotUnderPressured,obj.PressureBalanced,obj.Choked,obj.Flared,obj.Valid);
%             end
        end
        
    end
    methods (Static)
        
        function [obj,report] ...
                = checkConduitTolerances(cI,cO,obj)
            % [Z0pass,UnderPressurePass,ChokePass,PressBalancePass,FragCheck,FlareCheck,valid] = checkConduitResult(cO,Zthresh,Mthresh,Pthresh)
            % Check conduit output for validity
            % IN:   cO      = conduit model output struct
            %       Ztol = maximum conduit "top depth" (m). 
            %                   Default = 0.2*r0 (0.1 conduit diameter)
            %       Mtol = Check Mach number greater than minimum thresh = 1 - Mtol.
            %                   Default Mtol = 0.05;
            %       Ptol = Check abs((P_vent - P_surface)/P_surface)< Pthresh. 
            %                 ie check if vent and surface pressures differ by less
            %                 than this fraction.
            %                   Default = 0.002;
            %
            % PASS CONDITIONS ((ALWAYS REQ'D):
            %   Z0pass              1) TRUE if upper-most z value is less than Zthresh
            %                       ie Conduit run reaches surface. Typical Zthresh is 
            %                       within a conduit diameter or so
            %
            %   UnderPressurePass  2) TRUE if (P_v + rho*u^2/2)/P > (1- Pthresh) : Conduit is not
            %                       underpressured (ie likely to collapse).
            %
            %  AND EITHER OF
            % ChokePass             3) TRUE if Mach # < Mthresh with P_v/P_a > (1-Pthresh) :
            %                        Conduit is choked and either pressure balanced or overpressured
            % PressBalancePass    4) TRUE if Mach # < Mthresh and |dP/P| < (1 +/- Pthresh)
            %                       Conduit is not choked but is pressure balanced
            %
            % Additional diagnostics:
            % FragCheck    5) Did fragmentation occur?
            % FlareCheck   6) Does conduit flare? ie vent radius > initial radius
            %
            % valid = true/false for successful run based on the above conditions

            % if nargin<5 || isempty(verbose)
            %     verbose = false;
            % end
            % if nargin<4 || isempty(Ptol)
            %     Ptol = ...
            % end




            % Was phi_crit reached?
            phiCritReached = cO.porosity(end) >= cO.Par.phi_frag;
            % Did fragmentation occur?
            FragCheck   = logical(cO.Outcome.frag);

            % ---- Minimum depth tolerance ---- 
            zPass   = cO.Z(end)./cI.conduit_radius < cI.zFailTol; % Reaches surface

            % ---- Underpressure tolarance ---- 
                % 1 : total flow pressure (fluid + dynamic) ~ exit static pressure
            % UnderPressurePass  = ( (cO.pm(end)+cO.U(end).^2*cO.rho_magma(end)/2)/cO.Outcome.pf ) > (1-Pthresh); % Not underpressured

                % 2 : fluid pressure ~ exit static pressure 
            underPressurePass  =  ( 1 - (cO.pm(end))/cO.Par.pf )  < cI.pFailTol; % Not underpressured

            % ---- Pressure balance tolerance  ---- 
            pressBalancePass   = and( (cO.pm(end)/cO.Par.pf) < (1+cI.pFailTol), underPressurePass );              % Pressure balanced

            % ---- Mach number tolerance for choked case  ---- 
            % ChokePass   = and( cO.M(end)>Mthresh,cO.M(end)<1.1 );                           % Choked - upper limit a bit arbitrary here as it isn't really used
            chokePass   = abs( 1-cO.M(end)) < cI.mFailTol;

            % Fcheck  = cO.a(end) > cO.a(1);
            flareCheck  = abs((cO.a(end) - cO.a(1))./cO.a(1)) > .005;

            % Non - fragmented case corresponds to and(Zpass,UnderPressurePass) && PressBalancePass
            %   --> unless we allow flow to not reach surface...

            % Fragmented case
            valid   = and(zPass,underPressurePass) && or(pressBalancePass,chokePass);

            report = sprintf('R= %.5f, Z ~ 0: %i, Frag: %i, UnderP.: %i, P bal.: %i, Mach#: %i, Flare: %i, Valid: %i\n',...
                    cO.Par.a,zPass,FragCheck,underPressurePass,pressBalancePass,chokePass,flareCheck,valid);

            obj.ZFailTol        = cI.zFailTol;
            obj.MFailTol        = cI.mFailTol;
            obj.PFailTol        = cI.pFailTol;
            obj.DepthFlag       = zPass;
            obj.NotUnderPressured   = underPressurePass;
            obj.PressureBalanced    = pressBalancePass;
            obj.Valid           = valid;
            obj.Frag            = FragCheck;
            obj.PhiCritReached  = phiCritReached;
            obj.Choked           = chokePass;
            obj.Flared           = flareCheck;

%             outcome.Outcome = interpretConduitResult(outcome,verbose);

            % passThresh.underPressure = 

        end
        
        function name = interpretConduitResult(obj,verbose)
            % obj = struct of flags from checkConduitTolerances function
            % verbose = T/F to print outcome messages and additional error codes.
            %
            % OUT:
            %   code        = outcome code (int)
            %   outcome     = outcome name (char)
            %   msg         = detailed description (char)

            % if nargin<2 || isempty(pF) % Assume default thresholds
            %     cI = getConduitSource;
            %     Par.zFailTol = cI.zFailTol;
            %     Par.pFailTol = cI.pFailTol;
            %     Par.mFailTol = cI.mFailTol;
            % end

            if nargin<2 || isempty(verbose)
                verbose = false;
            end


            if obj.Valid && ~obj.Frag % Valid effusive
                if ~obj.PhiCritReached && obj.Valid  % Effusive
                    name = "validEffusive";

                elseif obj.PhiCritReached && obj.Valid % Effusive - ambiguous
                    name = "validEffusiveAmbiguous";
                else % Unmapped
                    name = "unmappedValidEffusive";
                end

            elseif ~obj.Valid && ~obj.Frag % Invalid effusive
                if obj.NotUnderPressured && ~obj.PressureBalanced && ~obj.Valid % Stalled/intrusive
                    name = "invalidEffusive";

                elseif obj.PressureBalanced % Effusive non-physical
                    name = "invalidIntrusive";
            %         code    = -2;
                else
                    name = "unamppedInvalidEffusive";
            %         code = UM.invalidEffusive;

                end

            elseif obj.Frag && obj.Valid % Valid explosive
                if obj.PressureBalanced % Pressure balanced jet
                    name = "validFragPressBalance";

                elseif ~obj.PressureBalanced && obj.Choked && ~obj.Flared % Overpressured, choked
                    name = "validExplosiveChoked";

                elseif ~obj.PressureBalanced && obj.Choked && obj.Flared % Overpressured, choked, flared
                    name = "validExplosiveFlaring";
                else
                    name = "unmappedValidExplosive";
                end

            elseif obj.Frag && ~obj.Valid % Invalid explosive
                if ~obj.DepthFlag && obj.PressureBalanced % Did not reach surface, pressure balanced
                    name = "invalidFragStalled";

                elseif ~obj.DepthFlag && obj.Choked && ~obj.Flared % Did not reach surface, M~1 limit reached immediately and without flaring (probably overpressure)
                    name = "invalidFragChoked";

                elseif ~obj.DepthFlag && obj.Choked && obj.Flared % Choked and flaring, M~1 limit reached (probably overpressure)
                    name = "invalidFragFlaring";

                elseif obj.DepthFlag && ~obj.Choked && ~obj.PressureBalanced % Reached surface, overpressured but did not reach M~1
                    name = "invalidFragNoChoke";
                else
                    name = "unmappedInvalidExplosive";
                end

            end

%             outcome = OutcomeCode(name,pF);

%             if verbose
%                 disp(obj)
%                 if ~isempty(obj.Exception)
%                     disp(obj.Exception)
%                     disp(cI)
%                 end
%             end

        end
        
        function T = getTable(obj)
            if nargin<1 || isempty(obj) % Assume default thresholds
                obj = ConduitOutcome.getDefaultTolerances;
            end
            
            names = ["nullResult"
                "validEffusive"
                "validEffusiveAmbiguous"
                "validFragPressBalance"
                "validExplosiveChoked"
                "validExplosiveFlaring"
                "invalidEffusive"
                "invalidIntrusive"
                "invalidFragStalled"
                "invalidFragChoked"
                "invalidFragFlaring"
                "invalidFragNoChoke"
                "unmappedValidEffusive"
                "unmappedValidExplosive"
                "unmappedInvalidEffusive"
                "unmappedInvalidExplosive"
                "failedIntegration"
                "failedPhysics"
                "failedUnmappedError"
                ];

%             codes = [0 1:5, -1:-1:-6, 21, 22, -21, -22, -99]';
            codes = [0 1:5, -1:-1:-6, 11, 12, -11, -12, -20, -21, -22]';


            label = ["Code 0 - no result"
                "Valid effusive"
                "Valid effusive - ambiguous"
                "Valid fragmenting - pressure balanced jet"
                "Valid explosive, choked at vent"
                "Valid explosive, choked and flaring"
                "Invalid effusive, non-physical"
                "Invalid stalled/intrusive"
                "Invalid fragmenting - stalled"
                "Invalid fragmenting - choked at Mach limit"
                "Invalid fragmenting - choked at Mach limit, flaring"
                "Invalid fragmenting - overpressured without choke"
                "UNMAPPED - valid effusive"
                "UNMAPPED - valid explosive"
                "UNMAPPED - invalid effusive"
                "UNMAPPED - invalid explosive"
                "FAILED, INTEGRATION TOLERANCE"
                "FAILED, PHYSICS ERROR"
                "FAILED, UNMAPPED ERROR"];

            msgs = { 'Code 0 - no result'
                'Valid effusive eruption.'
                'Fragmentation condition reached very near surface but no fragmentation step was performed. Valid to within tolerance.'
                sprintf('Fragmentation occurred, gas/pyroclast mixture reaches vent with Pm ~ Pf (Tol: %.2f).',obj.PFailTol)
                sprintf('Fragmentation occurred, choking occurs within %.1f vent radii of surface.',obj.ZFailTol)
                'Fragmentation occurred, gas/pyroclast choked with conduit flare towards vent.'
                'Effusive eruption was overpressured beyond tolerance at the vent.'
                'No fragmentation occurred, and pressure reached surface value far below the vent. Taken as invalid, but could be considered intrusive.'
                'Fragmentation occurred, but pressure reached surface value too far below the vent and integration stopped. Taken as invalid.'
                'Did not reach surface, M>0.999 limit reached immediately and without flaring (probably overpressured)'
                'Did not reach surface, M>0.999 limit reached after flaring (probably overpressure)'
                'Reached surface, overpressured without M~1 - insufficient decompression'
                'Valid, unfragmented conduit outcome - not previously identified - see ConduitOutcome.printTable.'
                'Valid, fragmenting conduit outcome - not previously identified - see ConduitOutcome.printTable.'
                'Invalid, unfragmented conduit outcome - not previously identified - see ConduitOutcome.printTable.'
                'Invalid, fragmenting conduit outcome - not previously identified - see ConduitOutcome.printTable.'
                'Simulation failed to meet integration tolerances.'
                'Simulation failed with known physics exception.'
                'Simulation failed with unmapped exception.'
                };

            T = table(codes, label, msgs, 'VariableNames', {'Code', 'Label', 'Message'}, 'RowNames', names);

        end
        
        function printTable(Tols)
            if nargin<1 || isempty(Tols) % Assume default thresholds
                Tols = getDefaultTolerances;
            end
            T = ConduitOutcome.getTable;
            disp(T)
        end
        
        function obj = getDefaultTolerances(obj)
            cI = getConduitSource;
            obj.ZFailTol = cI.zFailTol;
            obj.PFailTol = cI.pFailTol;
            obj.MFailTol = cI.mFailTol;
        end
        
        function cO = getErrorOutcomeFields
            % Populate NaN's in key fields to enable outcome plotting
            cO.Z        = NaN;
            cO.pm       = NaN;
            cO.M        = NaN;
            cO.porosity = NaN;
            cO.a        = NaN;
            cO.Par.pf   = NaN;
            cO.Par.phi_frag = NaN;
        end
    end
end

% function Tols = getDefaultTolerances
%     cI = getConduitSource;
%     Tols.zFailTol = cI.zFailTol;
%     Tols.pFailTol = cI.pFailTol;
%     Tols.mFailTol = cI.mFailTol;
% end

