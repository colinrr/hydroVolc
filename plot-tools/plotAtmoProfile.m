function [ax,p1] = plotAtmoProfile_pap(atmo,axi,col,plotWind)

if ischar(atmo)
    load(atmo)
    atmo = atmprofile;
end

if nargin<2
    axi = [];
end
if nargin<3
    col = [];
end
if nargin<4
    plotWind = true;
end

if plotWind
    nc = 4;
else
    nc = 3;
end

if isstruct(atmo)
    atmo = struct2table(rmfield(atmo,{'lat','lon','time','Units','years'}));
end

nr = 1;
% nc = 3;
dx = 0.01;
ppads = [0.1 0.03 0.1 0.07];
lw = 1.5;
fs = 9;

if isempty(axi)
    plotnew = true;
    if isempty(col)
        col = [ 0    0.4470    0.7410];
    end
else
    assert(length(axi)==(nr*nc))
    plotnew = false;
    ax = axi;
    if isempty(col)
        col = [0.8500    0.3250    0.0980];
    end
end

Z = atmo{:,'Altitude'}/1e3;
T = atmo{:,'Temperature'};
P = atmo{:,'Pressure'};

% %     V = nan(size(U));
% %     wl = {'Abs. val.'};
% else
%     V = atmo{:,'ZonalWindSpeed'};
%     U = atmo{:,'MeridionalWindSpeed'};
%     wl = {'U','V'};
% end
RH = atmo{:,'relativeHumidity'};

htropo = findTPheight(Z,T);

if plotnew
    figure('position',[50 50 1000 800])
end

for si = 1:(nr*nc)

    if plotnew
        ax(si) = tightSubplot(nr,nc,si,dx,[],ppads);
    end
    
    switch si
        case 1
            p1 = plot(ax(si),T,Z,'LineWidth',lw,'Color',col);
            xl = xlim(ax(si));
            text(ax(si),diff(xl)*0.85+xl(1),htropo+.2,'Tropopause',...
                'HorizontalAlignment','right','VerticalAlignment','bottom',...
                'FontSize',fs,'Color',col)
            if plotnew
                ylabel('Height (km a.s.l.)','Interpreter','Latex')
                xlabel('Temperature (K)','Interpreter','Latex')
            end

        case 2
%             ax(2) = tightSubplot(nr,nc,2,dx,[],ppads);
            plot(ax(si),P/1e3,Z,'LineWidth',lw,'Color',col)
            if plotnew
                xlabel('Pressure (kPa)','Interpreter','Latex')
            end
    % ax(3) = tightSubplot(nr,nc,3,dx,[],ppads);
    % plot(U,Z,'LineWidth',2)
    % hold on
    % plot(V,Z,'LineWidth',2)
    % xlabel('Wind (m/s)')
        case 3
%             ax(3) = tightSubplot(nr,nc,3,dx,[],ppads);
            plot(ax(si),RH,Z,'LineWidth',lw,'Color',col)
            if plotnew
                xlabel('Rel. Humidity (\%)','Interpreter','Latex')
            end
            
        case 4
            uv = [];
            wl = {};
            if ismember('MeridionalWindSpeed',atmo.Properties.VariableNames)
                uv = [uv atmo{:,'MeridionalWindSpeed'}];
                wl = [wl {'U'}];
            end
            if ismember('ZonalWindSpeed',atmo.Properties.VariableNames)
                uv = [uv atmo{:,'ZonalWindSpeed'}];
                wl = [wl {'V'}];
            end
            if ismember('windAbs',atmo.Properties.VariableNames)
                wabs = atmo{:,'windAbs'};
                wl = [wl {'Abs.'}];
            end
            if ~isempty(uv)
                plot(uv,Z,'LineWidth',2)
            end
                
            
            hold on
            plot(wabs,Z,'--k','LineWidth',2)
            xlabel('Wind (m/s)')
            legend([wl 'Tropopause'])
            
    % ax(5) = tightSubplot(nr,nc,5,dx,[],ppads);

    set(ax(2:nc),'YTickLabel',[])
    linkaxes(ax,'y')
    axis(ax,'tight')

    for ii=1:nc
        hold(ax(ii),'on')
        if ii==1
            tl = plot(ax(ii),xlim(ax(ii)),htropo*[1 1],'--','LineWidth',lw,'Color',col);
        else
            plot(ax(ii),xlim(ax(ii)),htropo*[1 1],'--','LineWidth',lw,'Color',col)
        end
    end
    % legend(ax(3),wl)
    % legend(tl,'Tropopause')
    xl = xlim(ax(1));
    set(ax,'FontSize',fs)

    end
    if plotnew
        for ai = 1:length(ax)
            hold(ax(ai),'on'); 
            grid(ax(ai),'on'); 
            set(ax(ai),'FontSize',fs)
%             xlabel(ax(ai),xl{ai},'Interpreter','Latex')
        end
        linkaxes(ax,'y')
        set(ax(2:end),'YTickLabel',{})
%         ylim(ax,wO.z([1 end]))
    end
end


% function plotAtmoProfile(atmo)
% 
% if ischar(atmo)
%     load(atmo)
%     atmo = atmprofile;
% end
% 
% nr = 1;
% nc = 4;
% dx = 0.01;
% ppads = [0.07 0.05 0.1 0.05];
% lw = 2;
% fs = 12;
% 
% Z = atmo{:,'Altitude'}/1e3;
% T = atmo{:,'Temperature'};
% P = atmo{:,'Pressure'};
% if ismember('Wind_abs',atmo.Properties.VariableNames)
%     U = atmo{:,'Wind_abs'};
%     V = nan(size(U));
%     wl = {'Abs. val.'};
% else
%     U = atmo{:,'Zonalwindspeed'};
%     V = atmo{:,'Meridionalwindspeed'};
%     wl = {'U','V'};
% end
% RH = atmo{:,'Relativehumidity'};
% 
% htropo = findTPheight(Z,T);
% 
% figure('position',[50 50 1000 800])
% ax(1) = tightSubplot(nr,nc,1,dx,[],ppads);
% plot(T,Z,'LineWidth',2)
% ylabel('Height (km a.s.l.)')
% xlabel('T (K)')
% 
% ax(2) = tightSubplot(nr,nc,2,dx,[],ppads);
% plot(P,Z,'LineWidth',2)
% xlabel('P (Pa)')
% 
% ax(3) = tightSubplot(nr,nc,3,dx,[],ppads);
% plot(U,Z,'LineWidth',2)
% hold on
% plot(V,Z,'LineWidth',2)
% xlabel('Wind (m/s)')
% 
% ax(4) = tightSubplot(nr,nc,4,dx,[],ppads);
% plot(RH,Z,'LineWidth',2)
% xlabel('Rel. Hum.')
% 
% % ax(5) = tightSubplot(nr,nc,5,dx,[],ppads);
% 
% set(ax(2:nc),'YTickLabel',[])
% linkaxes(ax,'y')
% axis(ax,'tight')
% 
% for ii=1:nc
%     hold(ax(ii),'on')
%     if ii==1
%         tl = plot(ax(ii),xlim(ax(ii)),htropo*[1 1],'--k','LineWidth',2);
%     else
%         plot(ax(ii),xlim(ax(ii)),htropo*[1 1],'--k','LineWidth',2)
%     end
% end
% legend(ax(3),wl)
% % legend(tl,'Tropopause')
% xl = xlim(ax(1));
% text(ax(1),diff(xl)*0.95+xl(1),htropo+.2,'Tropopause',...
%     'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fs)
% set(ax,'FontSize',fs)
% 
% 
% end
