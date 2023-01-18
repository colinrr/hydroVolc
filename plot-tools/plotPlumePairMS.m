function [axP,lp] = plotPlumePairMS(D,cols)
% plotCouplePlume(D)
% plots the detailed output of a single coupled conduit-plume run
%
% CRowell Mar 2021

    if nargin<2
        cols = get(0,'DefaultAxesColorOrder');
        cols = cols([2 1 3 4 5 6 7],:);
    end
    
    nrows = 1;
    ncols = 6;

    figure('position',[50 50 1200 600])

    % CONDUIT
    %   Pressure
    %   Velocity
    %   Decompression rate
    %   Gas volume fraction
    %   Nucleation
    
        % PLUME
    %   Density
    %   Velocity
    %   Radius
    %   Temperature
    %   Volume (mass?) fractions of particles, gas, liquid
    %   ...PSD?
    %   ...SO2?
    for ii=1:length(D)
        hm(ii) = D(ii).pO.hm;
    end
    [hm,pi] = max(hm);
%     if pi==2
%         d = d1; d1 = d2; d2 = d;
%     end
    
    if all(isfield(D,{'pI','pO'}))
        [axP,p1] = plotPlume(D(1).pI,D(1).pO,cols(1,:));
       
        for ii=2:length(D)
            p2 = plotPlume2(D(ii).pI,D(ii).pO,axP,cols(2,:));
        end
    end
    set(axP,'Ylim',[0 hm/1e3])
    lp = [p1 p2];
    
end


function [ax,p1] = plotPlume(pI,pO,color)



    % ylim
    yl = 6;
%     color = [0 0 .8];

    lw = 1.25; % Line Width  
    fs = 8; % Font size
    fscb = 8; %Legend font size
    
    dx = 0.02;
    ppads = [0.05 0.02 0.15 0.03];

%     figure
    nf = 6; % Number of subplots
    nc = 1;
    np = 0;
    
    ptp = ~isempty(pO.atmo.ztropo);
    
    % ===================================================
    % Density
    np = np + 1;
    ax(1) = tightSubplot(nc,nf,np,dx,[],ppads);
    hold on
    
    p1 = plot(pO.rho_B,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:));
    
    % Ambient
    yl = [0 max(pO.z)/1e3];
    plot(pO.atmo.rho,pO.z/1e3,'--k')
    xl = [0 5]; %xlim;
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw); end
    plot(xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))
    ax(np) = gca;
%     xlim([0 150])
    ylim(yl)
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Density, $\rho$ [kg/m$^3$]'},'Interpreter','Latex')
    ax(np).XTick = [0 5];
    ax(np).XMinorTick = 'on';
    ax(np).XMinorGrid = 'on';
    ax(np).MinorGridLineStyle = '-';
    %     ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    ylabel('Height, $z$ [km a.v.l.]','Interpreter','Latex')
    grid on
    
    text(4.75,pO.hb/1e3+.05,'$Z_{nbl}$','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fs,'Interpreter','Latex')
    text(4.75,pO.atmo.ztropo/1e3+.05,'Tropopause','HorizontalAlignment','right',...
        'VerticalAlignment','bottom','FontSize',fs,'FontName','Helvetica','Interpreter','tex')

    % ===================================================
    % Velocity
    np = np + 1;
    ax(2) = tightSubplot(nc,nf,np,dx,[],ppads);
    hold on
    
    plot(pO.u.*sin(pO.angle),pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
    
    % Ambient
%     plot(pO.rho_aB,pO.z/1e3,'--k')
    ax(np) = gca;
%     xlim([0 150])
    ylim([0 max(pO.z)/1e3])
    xl = [0 200]; %xlim;
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw); end
    plot(xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Velocity, $u$ [m/s]'},'Interpreter','Latex')
%     ax(np).XTick = [0:50:400];
%     ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on

    % ===================================================
    % Radius
    np = np + 1;
    ax(3) = tightSubplot(nc,nf,np,dx,[],ppads);
    hold on
    
    plot(pO.r/1e3,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
%     set(gca,'XScale','log') % YEH??

    % Ambient
%     plot(pO.rho_aB,pO.z/1e3,'--k')
    xl = [0 10];
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw); end
    plot(xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))

    ax(np) = gca;
    xlim(xl)
    ylim([0 max(pO.z)/1e3])
%     ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Radius, $a$ [km]'},'Interpreter','Latex')
%     ax(np).XTick = [0:50:400];
%     ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on    

    % ===================================================
    % Temperature
    np = np + 1;
    ax(4) = tightSubplot(nc,nf,np,dx,[],ppads);
    hold on
    
    plot(pO.theta,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
    
    % Ambient
    plot(pO.atmo.theta_a,pO.z/1e3,'--k')
    ax(np) = gca;
    xl = [0 1200]; %xlim;
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw); end
    plot(xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))

    %     xlim([0 150])
    ylim([0 max(pO.z)/1e3])
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Temperature, $T$ [K]'},'Interpreter','Latex')
%     ax(np).XTick = [0:50:400];
%     ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on   
    
    % ===================================================
    % Particle Mass fractions
    np = np + 1;
    ax(5) = tightSubplot(nc,nf,np,dx,[],ppads);
    hold on
    
    % Ambient vapour mixing ratio
%     plot(pO.atmo.w_a,pO.z/1e3,'--',...
%     'LineWidth',lw,'Color',color(1,:))
    % Saturation vapour mixing ratio 
%     plot(pO.atmo.w_s,pO.z/1e3,'--k',...
%     'LineWidth',lw)

    plot((pO.m_s)./pO.m,pO.z/1e3,...  % particle mass fraction
    'LineWidth',lw,'Color',color)
%     plot((pO.m_v)./pO.m,pO.z/1e3,...  % vapour mass fraction
%     '--','LineWidth',lw,'Color',color)
%     plot((pO.m_d)./pO.m,pO.z/1e3,...  % dry air mass fraction
%     'LineWidth',lw)    
%     plot((pO.m_l)./pO.m,pO.z/1e3,...  % liguid water mass fraction
%     'LineWidth',lw)    
    set(gca,'XScale','log')
    
    xl = [1e-2 1e0];
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw); end
    plot(xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))

    ax(np) = gca;
    xlim(xl)
    ylim([0 max(pO.z)/1e3])
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Particle Mass'; 'Fraction, $n_s$'},'Interpreter','Latex')
%     ax(np).XTick = [0:50:400];
%     ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on   
    
    % ===================================================
    % water mass fractions
    np = np + 1;
    ax(6) = tightSubplot(nc,nf,np,dx,[],ppads);
    hold on
    
%     xl = [1e-7 1e0];
%     pcolor(log10(pI.Rgsd),pO.z/1e3,pO.m_si)
%     plot((pO.m_s)./pO.m,pO.z/1e3,...  % particle mass fraction
%     'LineWidth',lw,'Color',color)
    plot((pO.m_v)./pO.m,pO.z/1e3,...  % vapour mass fraction
    'LineWidth',lw,'Color',color)
%     plot((pO.m_d)./pO.m,pO.z/1e3,...  % dry air mass fraction
%     'LineWidth',lw)    
    plot((pO.m_l)./pO.m,pO.z/1e3,...  % liguid water mass fraction
    '--','LineWidth',lw,'Color',color) 
    set(gca,'XScale','log')
    
    xl = [1e-4 1e0];
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw); end
    plot(xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))
    
%     if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw); end
%     plot(xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))
    pl(1) = plot(nan,nan,'LineWidth',lw,'Color',[0.3 0.3 0.3]);
    pl(2) = plot(nan,nan,'--','LineWidth',lw,'Color',[0.3 0.3 0.3]);
    [hh,icons,plots,txt] = legend(pl,{'Vapor, $n_v$','Liquid, $n_l$'},'location','northeast','FontSize',fscb,'Interpreter','Latex'); %,'FontName','Helvetica');
    icons(3).XData = [0.05 0.25];
    icons(5).XData = [0.05 0.25];

    ax(np) = gca;
    xlim(xl)
    ylim([0 max(pO.z)/1e3])
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Water Mass'; 'Fractions'},'Interpreter','Latex')
%     ax(np).XTick = [0:50:400];
%     ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on   

    for np = 2:(nf*nc)
        ax(np).YTickLabel = {};
%         ax(np).YTick = [0:20];
    end    
    
end


function p2 = plotPlume2(pI,pO,ax,color)


%     co = get(gca,'ColorOrder');
%     color = co(2,:);

    lw = 1.25; % Line Width  
    fs = get(gca,'FontSize');
    
    % ===================================================
    % Density
    p2 = plot(ax(1),pO.rho_B,pO.z/1e3,...
    'LineWidth',lw,'Color',color);
    xl = xlim(ax(1));
    plot(ax(1),xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))
    text(ax(1),4.75,pO.hb/1e3+.05,'$Z_{nbl}$','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fs,'Interpreter','Latex')

    % ===================================================
    % Velocity

    plot(ax(2),pO.u.*sin(pO.angle),pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
    xl = xlim(ax(2));
%     if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(ax(2),xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))
    
    % ===================================================
    % Radius
    plot(ax(3),pO.r/1e3,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
%     set(gca,'XScale','log') % YEH??

    % Ambient
%     plot(pO.rho_aB,pO.z/1e3,'--k')
    xl = xlim(ax(3)); %[1 1e4];
%     if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(ax(3),xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))
    
    % ===================================================
    % Temperature
    plot(ax(4),pO.theta,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
    
    % Ambient
%     plot(pO.atmo.theta_a,pO.z/1e3,'--k')
%     ax(np) = gca;
    xl = xlim(ax(4));
%     if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(ax(4),xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))

    % ===================================================
    % Mass or volume fractions?

    % Ambient vapour mixing ratio
%     plot(pO.atmo.w_a,pO.z/1e3,'--',...
%     'LineWidth',lw,'Color',color(1,:))
    % Saturation vapour mixing ratio 
%     plot(pO.atmo.w_s,pO.z/1e3,'--k',...
%     'LineWidth',lw)

    plot(ax(5),(pO.m_s)./pO.m,pO.z/1e3,...  % particle mass fraction
    'LineWidth',lw,'Color',color)
%     plot(ax(5),(pO.m_v)./pO.m,pO.z/1e3,...  % vapour mass fraction
%     '--','LineWidth',lw,'Color',color)
%     plot((pO.m_d)./pO.m,pO.z/1e3,...  % dry air mass fraction
%     'LineWidth',lw)    
%     plot((pO.m_l)./pO.m,pO.z/1e3,...  % liguid water mass fraction
%     'LineWidth',lw)    
%     set(gca,'XScale','log')
    
    xl = xlim(ax(5)); %[1e-4 1e0];
%     if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(ax(5),xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))
    
    % ===================================================
    % water mass fractions
    plot(ax(6),(pO.m_v)./pO.m,pO.z/1e3,...  % vapour mass fraction
    'LineWidth',lw,'Color',color)
%     plot((pO.m_d)./pO.m,pO.z/1e3,...  % dry air mass fraction
%     'LineWidth',lw)    
    plot(ax(6),(pO.m_l)./pO.m,pO.z/1e3,...  % liguid water mass fraction
    '--','LineWidth',lw,'Color',color)
    xl = xlim(ax(6));
    plot(ax(6),xl,pO.hb/1e3*[1 1],':','LineWidth',lw,'Color',color(1,:))

end
