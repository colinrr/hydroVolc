function [fig,ax,th,h] = plotConduitObjective(dat,xvar,ax,th,col)
%
% IN: dat   = struct or array of structs, each containing IO fields from
%           conduit model (cI,cO)
%     ax    = axes from a pre-existing call (to add new results to existing plot)
%     xvar  = optional string matching a field of the input struct cI,
%           which will be used as the x axis variable (i.e. for parameter
%           sweeps). Otherwise, x axis is run index. 
%            - xvar won't currently work for non-scalar input variables (eg
%            composition)
%            - multiple calls with overlapping xvars could get messy. Use
%            col input to mitigate
%     col   = specify plot color
%
% OUT: ax = axes handles
%      th = threshold line handles
%      h  = scatter handles
%
% Temporary function version prior to building into generalized sweep
% function?
% To-do: include plots of stoppage/transitions criteria variables
%   z already in
%   phi
%   BND?
%   ...

assert(isstruct(dat),'Input D must be struct with fields cI,cO')
Nruns = length(dat);

if nargin<3 || isempty(ax)
    use_new_axes = true;
%     ax = [];
elseif isa(ax,'matlab.graphics.axis.Axes')
    use_new_axes = false;
end
% if nargin<3 || isempty(th)
% end
if nargin<2 || isempty(xvar)
    xvar = [];
    use_x_var = false;
    xlab = 'Run No.';
elseif isfield(dat(1).cI,xvar)
    use_x_var = true;
    xvals = zeros(Nruns,1);
    xlab = xvar;
end
if nargin<5
    col = [0 0 0]; %get(0,'DefaultAxesColorOrder');
%     col = col(1,:);
end

% Plots setup
% -- With top colorbar
% figpos = [100 100 950 1000];
% ppads = [0.15 0.05 0.06 0.17];
% cbpos = [0.85 0.025];
% cbposlabel = 'southoutside';

% -- With side colorbar
figpos = [100 100 1050 900];
ppads = [0.14 0.2 0.06 0.02];
cbpos = [0.82 0.025];
cbposlabel = 'eastoutside';
dx = [];
dy = 0.013;

msz = 100;
msym = 's';
lw = 1.0;
fs = 15;
fscb = 11;
lstyle = '--';

idx1 = 0;
if use_new_axes
    fig = figure('position',figpos);
    nr = 7;
    nc = 1;
    
    clear ax
    za = 1; ax(za) = tightSubplot(nr,nc,za,dx,dy,ppads); % Z
    pa = 2; ax(pa) = tightSubplot(nr,nc,pa,dx,dy,ppads); % P
    ma = 3; ax(ma) = tightSubplot(nr,nc,ma,dx,dy,ppads); % M
    pha = 4; ax(pha) = tightSubplot(nr,nc,pha,dx,dy,ppads); % phi
    fla = 5; ax(fla) = tightSubplot(nr,nc,fla,dx,dy,ppads); % Flaring
    va = 6; ax(va) = tightSubplot(nr,nc,va,dx,dy,ppads); % Valid?: and(Z,~UP), or(choke,Pbal)
    fa = 7; ax(fa) = tightSubplot(nr,nc,fa,dx,dy,ppads); % Flags: frag, choke, flare
%     oa = 8; ax(oa) = subplot(nr,nc,oa); % Interpreted outcome
    hold(ax,'on')
    
    if ~use_x_var
        xvals = 1:Nruns;
    end
else
        
    if ~use_x_var % Get the last x index plotted
        ch = get(ax(1),'Children');
        xl = get(ax(1),'XLim');
        for hi = 1:length(ch)
            idx1 = max([idx1 max(ch(hi).XData(:))]);
        end
        xvals = idx1 + (1:Nruns);
    end
end


% build vecs
Zvec = zeros(Nruns,1);
Pvec = Zvec;
Mvec = Zvec;
PHvec = Zvec;
Flvec = Zvec;
% Zpass = false(Nruns,1);
% UPpass = false(Nruns,1);
% PBpass = false(Nruns,1);

failed             = false(Nruns,1);
z_and_not_up       = false(Nruns,1);
balance_or_choke   = false(Nruns,1);
frag                = false(Nruns,1);
choke               = false(Nruns,1);
flare               = false(Nruns,1);
valid               = false(Nruns,1);
outcomeCode         = zeros(Nruns,1);

zThresh = zeros(Nruns,1);
pThresh = zeros(Nruns,2);
mThresh = zeros(Nruns,2);
phThresh = ones(Nruns,1);
for nn = 1:Nruns
    if use_x_var
        xvals(nn) = dat(nn).cI.(xvar);
    end
    try
        Outcome = dat(nn).cO.Outcome;
    catch ME
        pause(0.5)
    end
    
    zThresh(nn) = Outcome.ZFailTol; %.*dat(nn).cI.conduit_radius;
%     Pthresh(nn,:) = [1-dat(nn).cI.pFailTol-(dat(nn).cO.U(end).^2*dat(nn).cO.rho_magma(end)/2)/dat(nn).cO.Par.pf 1+dat(nn).cI.pFailTol]; %./dat(nn).cI.pFailthresh;
    pThresh(nn,:) = 1 + [-1 1].*dat(nn).cI.pFailTol;
    mThresh(nn,:)   = 1+Outcome.MFailTol.*[-1 1]; % or just 1
%     [Zpass,UPpass,choke(nn),PressBalPass,frag(nn),flare(nn),valid(nn),~] 
%     outcome     = checkConduitTolerances(dat(nn).cI,dat(nn).cO);
    failed(nn)  = Outcome.Failed;
    choke(nn)   = Outcome.Choked;
    frag(nn)    = Outcome.Frag;
    flare(nn)   = Outcome.Flared;
    valid(nn)   = Outcome.Valid;
    outcomeCode(nn) = Outcome.Code;

    z_and_not_up(nn) = and(Outcome.DepthFlag,Outcome.NotUnderPressured);
    balance_or_choke(nn) = or(choke(nn),Outcome.PressureBalanced);
    Zvec(nn) = abs(dat(nn).cO.Z(end))./dat(nn).cI.conduit_radius; %Zthresh(nn);
    Pvec(nn) = dat(nn).cO.pm(end)./dat(nn).cO.Par.pf;
    Mvec(nn) = dat(nn).cO.M(end);
    PHvec(nn) = dat(nn).cO.porosity(end)./dat(nn).cO.Par.phi_frag;
    Flvec(nn) = max(dat(nn).cO.a)./dat(nn).cI.conduit_radius;
end





% Potential stuff to plot continuous threshold lines - may not play well
% with specified xvar?
%     lh2 = findobj(ax(2),'LineStyle',lstyle);
%     lh3 = findobj(ax(ma),'LineStyle',lstyle);
% if ~use_new_axes
%     lh1 = findobj(ax(1),'LineStyle',lstyle);
%     xdat = [lh1.XData xvals];
%     ydat = [lh1.YData Zthresh];
%     
% else
%     xdat = xvals;
%     ydat = Zthresh;
% end

%     failmap = [0.8 0.2 0.2; 0 0.6 1]; cax = [0 1];
[failmap,cax,cticks,cticklabels] = outcomeColorMap;
%     colorval = outcomeCode;
[~,colorval] = ismember(outcomeCode,sort(ConduitOutcome.getTable.Code));

checkmap = [0 0 0; 0 0.6 1];
validmap = [0.3 0 0; checkmap];


th(1) = plot(ax(za),xvals,zThresh,lstyle,'Color',col,'Linewidth',lw);
h(1)  = scatter(ax(za),xvals,Zvec,msz,colorval,'filled','MarkerEdgeColor',col,'LineWidth',lw,'Marker',msym);

th(2:3) = plot(ax(pa),xvals,pThresh,lstyle,'Color',col,'Linewidth',lw);
h(2)  = scatter(ax(pa),xvals,Pvec,msz,colorval,'filled','MarkerEdgeColor',col,'LineWidth',lw,'Marker',msym);

th(4:5) = plot(ax(ma),xvals,mThresh,lstyle,'Color',col,'Linewidth',lw);
h(3)  = scatter(ax(ma),xvals,Mvec,msz,colorval,'filled','MarkerEdgeColor',col,'LineWidth',lw,'Marker',msym);

th(6) = plot(ax(pha),xvals,phThresh,lstyle,'Color',col,'Linewidth',lw);
h(4)  = scatter(ax(pha),xvals,PHvec,msz,colorval,'filled','MarkerEdgeColor',col,'LineWidth',lw,'Marker',msym);

h(5)  = scatter(ax(fla),xvals,Flvec,msz,colorval,'filled','MarkerEdgeColor',col,'LineWidth',lw,'Marker',msym);

% xl = xlim(ax(va));
validIm = double([balance_or_choke z_and_not_up])';
validIm(:,failed) = -1;
h(6) = imagesc(ax(va),xvals,1:2,validIm);

xl = xlim(ax(va));
h(7) = imagesc(ax(fa),xvals,1:3,[flare frag choke]');

% h(8) = imagesc(ax(fa),xvals,1,outcomeCode);

all_ax = [za pa ma pha fla va fa];
if use_new_axes
    
%     colormap(ax(za),failmap)
%     colormap(ax(pa),failmap)
%     colormap(ax(ma),failmap)
%     colormap(ax(pha),failmap)

    for aa=[za pa ma pha fla]; caxis(ax(aa),cax); colormap(ax(aa),failmap); end
    colormap(ax(va),validmap)
    caxis(ax(va),[-1.5 1.5])
    colormap(ax(fa),checkmap)
    caxis(ax(fa),[-0.5 1.5])
    ylabs = {'$|Z_{min}|/a$', '$P_m/P_f$', '$M$', '$\phi/\phi_{frag}$', '$a_{max}/a_0$', 'Valid?', 'Flags'};
    for ai = 1:length(all_ax); ylabel(ax(all_ax(ai)),ylabs{ai},'Interpreter','latex'); end
    %     caxis(ax(za),[0 1])
%     caxis(ax(pa),[0 1])
%     caxis(ax(ma),[0 1])
%     caxis(ax(pha),[0 1])
%     caxis(ax(pha),[0 1])
%     caxis(ax(va),[0 1])
%     caxis(ax(fa),[0 1])
    axis(ax(va),'tight')
    axis(ax(fa),'tight')
%     set(ax([1 2 3]),'Yscale','log')
%     ylabel(ax(za),'$|Z_{min}|/a$','Interpreter','latex')
%     ylabel(ax(pa),'$P_m/P_f$','Interpreter','latex')
%     ylabel(ax(ma),'$M$','Interpreter','latex')
%     ylabel(ax(pha),'$\phi/\phi_{frag}$','Interpreter','latex')
%     ylabel(ax(fla),'$a_{max}/a_0$','Interpreter','latex')
%     ylabel(ax(va),'Valid?','Interpreter','latex')
%     ylabel(ax(fa),'Flags','Interpreter','latex')
    xlabel(ax(fa),replace(xlab,'_',' '),'Interpreter','latex')

    set(ax(va),'YTick',(1:2),'YTickLabel',{'or(choke, P-bal)','and(Z, ~U-P)'})
    set(ax(fa),'YTick',(1:3),'YTickLabel',{'Flaring','Fragmented','Choked'})
    
    grid(ax(1:5),'on')
    xlim(ax,[min(xvals)-0.5 max(xvals)+0.5])
    set(ax,'FontSize',fs)
    linkaxes(ax,'x')
    
    % Top colorbar
%     axpos = get(ax(fla),'Position');
%     cb = colorbar(ax(fla),'location',cbposlabel);
%     cb.Position = [axpos(1) cbpos(1) axpos(3) cbpos(2)];
%     cb.Ruler.TickLabelRotation=50;
    
    % Side colorbar
    axpos1 = get(ax(za),'Position');
    axpos2 = get(ax(fla),'Position');
    cb = colorbar(ax(fla),'location',cbposlabel);
    cb.Position = [cbpos(1) axpos2(2) cbpos(2) sum([axpos1([2 4]) -axpos2(2)])];
%     cb.Ruler.TickLabelRotation=-10;
    
    cb.Ticks = cticks;
    cb.TickLabels = cticklabels;
    cb.FontSize = fscb;
    
    % Small colorbars
    axpos1 = get(ax(va),'Position');
    cb = colorbar(ax(va),'location',cbposlabel);
    cb.Position = [cbpos(1) axpos1(2) cbpos(2) axpos1(4)];
    cb.Ticks = [-1 0 1];
    cb.TickLabels = {'Failed','Invalid','Valid'};
    cb.FontSize = fscb;
    
    axpos2 = get(ax(fa),'Position');
    cb = colorbar(ax(fa),'location',cbposlabel);
    cb.Position = [cbpos(1) axpos2(2) cbpos(2) axpos2(4)];
    cb.Ticks = [0 1];
    cb.TickLabels = {'False','True'};
    cb.FontSize = fscb;


else
%     xlim(ax,[xl(1) max(xvals)+0.5])
%     xlim(ax(va),[min(xvals)-0.5 max(vals)+0.5])
end

set(ax,'XTickMode','auto','XLimMode','auto')
if length(get(ax(fa),'XTick')) > (idx1+Nruns) && ~use_x_var
    set(ax,'XTick',1:Nruns)
end
set(ax(1:end-1),'XTick',{})
    
end

function [cmap,cax,cticks,clabels] = outcomeColorMap
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

% function th = plotThreshLine(ax,new_x, new_y)
%  To re-plot a continuous threshold in subsequent plots
%     
%     lh = findobj(ax,'LineStyle',lstyle);
%     
%     if ~isempty(lh)
%         for li = 1:length(lh)
%             
%         end
%     else
%     end
% 
% end