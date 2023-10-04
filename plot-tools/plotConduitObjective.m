function [ax,th,h] = plotConduitObjective(dat,xvar,ax,th,col)
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
% function

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

lstyle = '--';

idx1 = 0;
if use_new_axes
    figure('position',[100 100 950 700])
    nr = 5;
    nc = 1;
    ppads = [];
    
    clear ax
    ax(1) = subplot(nr,nc,1); % Z
    ax(2) = subplot(nr,nc,2); % P
    ax(3) = subplot(nr,nc,3); % M
    ax(4) = subplot(nr,nc,4); % Valid?: and(Z,~UP), or(choke,Pbal)
    ax(5) = subplot(nr,nc,5); % Flags: frag, choke, flare
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
% Zpass = false(Nruns,1);
% UPpass = false(Nruns,1);
% PBpass = false(Nruns,1);

Z_and_not_UP       = false(Nruns,1);
balance_or_choke   = false(Nruns,1);
frag = false(Nruns,1);
choke = false(Nruns,1);
flare = false(Nruns,1);
valid = false(Nruns,1);

Zthresh = zeros(Nruns,1);
Pthresh = zeros(Nruns,2);
Mthresh = zeros(Nruns,2);

for nn = 1:Nruns
    if use_x_var
        xvals(nn) = dat(nn).cI.(xvar);
    end
    
    Zthresh(nn) = dat(nn).cI.ZfailScale; %.*dat(nn).cI.conduit_radius;
    Pthresh(nn,:) = [1-dat(nn).cI.Pfailthresh-(dat(nn).cO.U(end).^2*dat(nn).cO.rho_magma(end)/2)/dat(nn).cO.Par.pf 1+dat(nn).cI.Pfailthresh]; %./dat(nn).cI.Pfailthresh;
    Mthresh(nn,:)   = [dat(nn).cI.Mfailthresh 2-dat(nn).cI.Mfailthresh]; % or just 1
    [Zpass,UPpass,choke(nn),PressBalPass,frag(nn),flare(nn),valid(nn)] = checkConduitResult(dat(nn).cO,Zthresh(nn).*dat(nn).cI.conduit_radius,dat(nn).cI.Mfailthresh,dat(nn).cI.Pfailthresh);

    Z_and_not_UP(nn) = and(Zpass,UPpass);
    balance_or_choke(nn) = or(choke(nn),PressBalPass);
    Zvec(nn) = abs(dat(nn).cO.Z(end))./dat(nn).cI.conduit_radius; %Zthresh(nn);
    Pvec(nn) = dat(nn).cO.pm(end)./dat(nn).cO.Par.pf;
    Mvec(nn) = dat(nn).cO.M(end);
end

% Plots
msz = 40;
lw = 1.0;
fs = 14;

% Potential stuff to plot continuous threshold lines - may not play well
% with specified xvar?
%     lh2 = findobj(ax(2),'LineStyle',lstyle);
%     lh3 = findobj(ax(3),'LineStyle',lstyle);
% if ~use_new_axes
%     lh1 = findobj(ax(1),'LineStyle',lstyle);
%     xdat = [lh1.XData xvals];
%     ydat = [lh1.YData Zthresh];
%     
% else
%     xdat = xvals;
%     ydat = Zthresh;
% end

th(1) = plot(ax(1),xvals,Zthresh,lstyle,'Color',col,'Linewidth',lw);
h(1)  = scatter(ax(1),xvals,Zvec,msz,valid,'filled','MarkerEdgeColor',col,'LineWidth',lw);

th(2:3) = plot(ax(2),xvals,Pthresh,lstyle,'Color',col,'Linewidth',lw);
h(2)  = scatter(ax(2),xvals,Pvec,msz,valid,'filled','MarkerEdgeColor',col,'LineWidth',lw);

th(4:5) = plot(ax(3),xvals,Mthresh,lstyle,'Color',col,'Linewidth',lw);
h(3)  = scatter(ax(3),xvals,Mvec,msz,valid,'filled','MarkerEdgeColor',col,'LineWidth',lw);

% xl = xlim(ax(5));
h(4) = imagesc(ax(4),xvals,1:2,[balance_or_choke Z_and_not_UP]');

xl = xlim(ax(5));
h(5) = imagesc(ax(5),xvals,1:3,[flare frag choke]');


if use_new_axes
    failmap = [0.8 0.2 0.2; 0 0.6 1];
    checkmap = [0 0 0; 0 0.6 1];
    colormap(ax(1),failmap)
    colormap(ax(2),failmap)
    colormap(ax(3),failmap)
    colormap(ax(4),checkmap)
    colormap(ax(5),checkmap)
    caxis(ax(1),[0 1])
    caxis(ax(2),[0 1])
    caxis(ax(3),[0 1])
    caxis(ax(4),[0 1])
    caxis(ax(5),[0 1])
    axis(ax(4),'tight')
    axis(ax(5),'tight')
    set(ax([1 2 3]),'Yscale','log')
    ylabel(ax(1),'$|Z_{min}|/a$','Interpreter','latex')
    ylabel(ax(2),'$P_m/P_f$','Interpreter','latex')
    ylabel(ax(3),'$M$','Interpreter','latex')
    ylabel(ax(4),'Valid?','Interpreter','latex')
    ylabel(ax(5),'Flags','Interpreter','latex')
    xlabel(ax(5),xlab,'Interpreter','latex')

    set(ax(4),'YTick',(1:2),'YTickLabel',{'or(choke, P-bal)','and(Z, ~U-P)'})
    set(ax(5),'YTick',(1:3),'YTickLabel',{'Flare','Frag','Choke'})
    
    grid(ax(1:3),'on')
    xlim(ax,[min(xvals)-0.5 max(xvals)+0.5])
    set(ax,'FontSize',fs)
    linkaxes(ax,'x')
else
    xlim(ax,[xl(1) max(xvals)+0.5])
%     xlim(ax(5),[min(xvals)-0.5 max(vals)+0.5])
end

set(ax,'XTickMode','auto')
if length(get(ax(5),'XTick')) > (idx1+Nruns) && ~use_x_var
    set(ax,'XTick',1:Nruns)
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