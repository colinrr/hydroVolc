function [cmap,cax,cticks,clabels,outcomeIndex,cbh] = outcomeColorMap(outcomeCodes, simplify, makeColorBar)
% [cmap,cax,cticks,clabels,outcomeIndex,cbh] = outcomeColorMap(simplify, makeLegend)
% OPTIONAL IN:
%   outcomeCodes: a set of outcome codes to map into the colormap. Makes
%           for easier plotting since outcome codes are not linear + monotonic      
%
%   simplify : boolean, if TRUE produces a simplify set of colours based on
%               on the hierarchy described below. FALSE produces a full set
%               of color for all established outcome codes.
%
%   makeColorBar : Enter as Boolean or any of the standard colorbar locations 
%               , eg 'northoutside'. Automatically generates colorbar 
%               on the current axes with appropriate orientation and labels.
%               Outputs handle.
%               Default location is 'eastoutside'.
%               Additional parameters such as font size, exact position,
%               etc must still be handled manually.
%
% OUT:
%   cmap            : the color map array
%   cax             : color axes limits appropriate to match the colormap
%   cticks          : colorbar tick positions
%   clabels         : colorbar tick labels
%   outcomeIndex    : color indices to use for plotting with cmap
%   cbh             : colorbar handle (if makeColorBar = true)
%
% C Rowell Jan 2023
%
% ---- Simple color hierarchy:-----
%  4 > 5 --> choked, overpressured explosive
%  3 --> pressure balance jet
%  2,1 --> effusive
%
% -2 --> 'intrusive'?
% Any other -ve - Invalid
% All Failed

% So 6 colors
%   - orange
%   - yellow
%   - blue

% - brown?
% - gray
% - black
% ---------------------------------

narginchk(0,3)

if nargin<1
    outcomeCodes = [];
end

if nargin<2 || isempty(simplify)
    simplify = false;
end
if nargin<3 || isempty(makeColorBar)
    makeColorBar = false;
end

assert(or(islogical(makeColorBar),ischar(makeColorBar)),...
    'Input "makeColorBar" must be of type logical or char.')

if simplify
    [cmap, cax, cticks, clabels, outcomeIndex] = simpleColorMap(outcomeCodes);
    
else
    if ~isempty(outcomeCodes)
        [~,outcomeIndex] = ismember(outcomeCodes,sort(ConduitOutcome.getTable.Code));
    end
    [cmap, cax, cticks, clabels] = fullColorMap;
end

if makeColorBar
    if islogical(makeColorBar)
        makeColorBar = 'eastoutside';
    end
    
    colormap(gca,cmap)
    cbh = colorbar(gca,'location',makeColorBar);
    caxis(cax);
    cbh.Ticks = cticks;
    cbh.TickLabels = clabels;
else
    cbh = [];
end

end

function [cmap, cax, cticks, clabels] = fullColorMap
    codes = ConduitOutcome.getTable.Code;
    [~,idx] = sort(ConduitOutcome.getTable.Code);
    for ii = 1:length(codes)
        clabels{ii} = sprintf('%i %s',codes(idx(ii)),ConduitOutcome.getTable.Properties.RowNames{idx(ii)});
    end
%     clabels = ConduitOutcome.getTable.Properties.RowNames(idx);
%     clabels = string(sort(codes));
    cticks = 1:length(codes);
    cax = [0.5 length(codes)+0.5];
    
    
    % Divides: [fail invalidUnmapped invalid 0(null) valid validUnmapped]
    codecut = [-20 -10 -1 1 10];
    
    % SETTING COLORMAPS FOR ALL OUTCOMES
    nullCol = [0.8 0.8 0.6];
%     validCols = [0 0 1; 0 0.8 0.8];
%     validCols1 = [0    0.5000    0.4000; 0.5000    0.7500    0.4000];
%     validUMcols = [0.0    0.9839    0.0805; 0.0    0.7993    0.3480];
%     validUMcols = [0.8333    0.5208    0.3317; 1.0000    0.7812    0.4975];
    validUMcols = [ 0.5250    0.6500    0.6500; 0.7625    0.8250    0.8250];
%     invalidCols = [1 1 0; 1 0.7 0.4];
%     invalidUMCols = [1 0 1; 0.8 0.3 0.8];
    invalidUMCols = [0.8 0.8 0.8; 0.6 0.6 0.6];
%     failCols = [1 0 0; 0.5 0 0];
    failCols = [0 0 0; 0.4 0.1 0.1];
    
    idxValid        = and( codes > 0, codes < codecut(5));
    idxValidUM      = codes >= codecut(5);
    idxInvalid      = and( codes >= codecut(2), codes <= codecut(3));
    idxInvalidUM    = and( codes > codecut(1), codes <= codecut(2));
    idxFail         = codes <= codecut(1);
    
    nValid = sum(idxValid);
    nValidUM = sum(idxValidUM);
    nFail  = sum(idxFail);
    nInvalid = sum(idxInvalid);
    nInvalidUM = sum(idxInvalidUM);
    
    
    cmap = [interpCols(validUMcols,nValidUM)
            winter(nValid)
%             interpCols(validCols,nValid)
            nullCol
            flipud(plasma(nInvalid))
%             interpCols(invalidCols,nInvalid)
            interpCols(invalidUMCols,nInvalidUM)
            interpCols(failCols,nFail)
            ];
    cmap = flipud(cmap);
        
%     validCMap = winter(nValid);
%     validUMcmap = 
%     invalidCMap = spring(nInvalid);
    
end

function [cmap, cax, cticks, clabels, outcomeIndex] = simpleColorMap(outcomeCodes)
% Simplified codes for quick visualization
% ---- Simple color hierarchy:-----
%  4 > 5 --> choked, overpressured explosive
%  3 --> pressure balance jet
%  2,1 --> effusive
%
% -2 --> 'intrusive'?
% Any other -ve - Invalid
% All Failed

% So 6 colors
%   - blue - explosive
%   - light blue - pressure balanced
%   - green - effusive

% - orange - intrusive
% - gray - other known/uknown invalid
% + STANDARD error codes
% ---------------------------------
narginchk(0,1)

if nargin<1
    outcomeCodes = [];
end


newCodes      = [-5:1:3];
oldCodeGroups = { -22, -21, -20, [-19:-3 0], -1, -2 , [1 2], 3, [4 5]};

clabels = {'Unmapped Error',
           'Failed Physics',
           'Failed Integration',
           'Other Invalid',
           'Invalid Effusive',
           'Intrusive',
           'valid Effusive',
           'valid Pressure Balanced'
           'valid Explosive',
           };

cax = [-5.5 3.5];
cticks = newCodes;
       
cmap = [0.6 0.1 0.1
        0.3 0.05 0.05
        0  0  0
        0.4 0.4 0.4
        0.8 0.8 0.5
        0.95 0.5 0.3
        0 0.8 0.4
        0 0.3 0.8
        0 0 0.6
    ];

outcomeIndex = zeros(size(outcomeCodes));

for ci = 1:length(newCodes)
    outcomeIndex(ismember(outcomeCodes,oldCodeGroups{ci})) = newCodes(ci);
    
end
end

function cmap = interpCols(cols,n)
    cmap = zeros(n,3);
    for ci = 1:3
%         if n==1
%             
%         else
            cmap(:,ci) = linspace(cols(1,ci),cols(2,ci),n)';
%         end
    end

end