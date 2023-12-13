function dat = conduitParameterSweep(cI_all,sweepParams,outDir,descriptor) % runCoarse,runFine) % OR conduitPressureSweep?
% Conduit sweep function, conduit V7
% INPUT:
%    cI = fixed conduit input params
%
%
%  ---> Shooting runs to find the range of viable conduit solutions for
%  "dry" runs, aka Hw=0
%       --> Conduit stop conditions: Pm~Pf, U~0, M>~.999
%       --> Check for valid conduit outcomes
%       --> Require Z(end)~0 AND either (a) choked or (b) pm(end) = pf
%           --> So watch out for runs where stop condition is M>=.999
%  ---> Test a range of vent pressures (~water depths) on the base case to
%   assess response, adjusting MER or radius (phi0 or Z0 as options?)
%
% PROCEDURE (should write as accessible functions to refine full look-up table):
%     (1) "runCoarse": Coarse brute force sweep with a coarse step size to find
%           approximate bounding locations for conduit radii.
%           --> Generate a coarse lookup table in R,MER,Pf
%           --> CURRENTLY SAVES ALL CONDUIT OUTPUT, FILE CAN BE VERY LARGE
%
%     (2) "getParams": Get summary values from coarse search, sufficient to
%           run fine search.
%           --> Could combine with (1) if you don't care to save coarse
%           output
%
%     (3) "runRefine": Use initial coarse table to:
%           (a) refine precision bounds on R for each MER,Pf 
%               --> (or go the other way and refine precision MER for given R,Pf)?
%
%   NOT YET FULLY IMPLEMENTED
%           (b) interpolate guesses for arbitrary R, MER, Pf from fine
%           search lookup table.
%           (c) run a fine search down to a tolerance from interpolations.
%               '-> This should be a fast-ish function to do a final,
%                   repeatable shoot (e.g. get MER for given input R, Pf)
%
%
% C Rowell, December 2023
% Uses Hajimirza conduit model, updated to V7



end

% Function to generate randomized model input
function [cI,pI,randPars] = populateParams(cI,SP)
% Populate input structs using sweep parameters in input struct
% randPars = vector of randomized variable results


    
    sn = {'cI'};
    ff{1} = fieldnames(SP.cI);
%     ff{2} = fieldnames(SP.pI);
    
    fo{1} = cI;
%     fo{2} = pI;
    
    randPars = struct();
    fungroups = cell(1,2); %Input function group
    
    for ss = 1:length(ff)
        if ~isempty(ff{ss})
            for ii = 1:length(ff{ss})
                fn = ff{ss}{ii};

    %             switch fn 
    %                 case 'logQ_rng'
    %                     fo{1}.Q = rand.*range(MC.logQ_rng) + min(MC.logQ_rng);
    %                     
    %                 case 'a_var'
    %                     a_peak = extrapVentRadius(cI.Q);
    %                     fo{1}.conduit_radius = a_peak.*(MC.a_var.*(2.*rand-1) + 1);
    %                     
    %                 otherwise

                assert(all(ismember({'dist','range'},fieldnames(SP.(sn{ss}).(fn)))),...
                    'Monte Carlo input: Variable "%s" is missing a required field',fn)
                dist = SP.(sn{ss}).(fn).dist;
                par  = SP.(sn{ss}).(fn).range;


                switch dist
                    case 'uniform'
                        fo{ss}.(fn) = range(par).*rand + min(par);

                    case 'discrete'
                        fo{ss}.(fn) = par(randi(length(par)));

                    case 'normal'
                        fo{ss}.(fn) = randn.*par(2) + par(1);

                    case 'function'
                        varfun.var = fn;
                        varfun.fun = par;
                        fungroups{ss} = [fungroups{ss} varfun];
                        fo{ss}.(fn) = [];
%                         fungroups{ss} = [fungroups{ss} {par}];
                        
                    otherwise
%                         if isa(par,'function_handle')
%                             fungroups{ss}.fun = par;
%                             fungroups{ss}.var = dist;
%                             tvar = eval(dist);
%                             fo{ss}.(fn) = range(tvar);
%                         else
                            warning('Distribution type "%s" not recognized, variable "%s" will not be randomized.',SP.(sn{ss}).(fn).dist,fn)
                            fo{ss}.(fn) = mean(par);
%                         end
                end        

                randPars.(fn) = fo{ss}.(fn);

            end
        end
    end

    % Special mass flux and conduit radii params
    % Mass flux
    if isfield(fo{1},'logQ')
        fo{1}.Q = 10^fo{1}.logQ;
        fo{1} = rmfield(fo{1},'logQ');
        
        randPars.Q = fo{1}.Q;
        randPars = rmfield(randPars,'logQ');
    end
    % Conduit radius lookup
    if isfield(fo{1},'a_var')
        a_peak = extrapVentRadius(fo{1}.Q);
        fo{1}.conduit_radius = a_peak.*fo{1}.a_var;
        fo{1} = rmfield(fo{1},'a_var');
        
        randPars.conduit_radius = fo{1}.conduit_radius;
        randPars = rmfield(randPars,'a_var');

    end
    
    cI = fo{1};
    pI = fo{2};
    
    % Run any input params that are a function of random params
    if ~isempty(fungroups{1})
        for jj=1:length(fungroups{1})
            cI.(fungroups{1}(jj).var) = eval(fungroups{1}(jj).fun);
%             cI.(fungroups{1}(jj).var) = varfun();
        end
    end
    if ~isempty(fungroups{2})
        for jj=1:length(fungroups{2})
            pI.(fungroups{2}(jj).var) = eval(fungroups{2}(jj).fun);
        end
    end
%     randPars = struct2table(randPars);
    


end