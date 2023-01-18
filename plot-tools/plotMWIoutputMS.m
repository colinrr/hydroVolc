% [ax,p1] = plotMWIoutput_pap(D,axi,color)

function [ax,p1] = plotMWIoutputMS(wO,cI,axi,color)

    % U, psi, T, rho, n_l, n_v, A
    nr = 1;
    nc = 6;
    dx = 0.02;
    dy = [];
    ppads = [0.06 0.02 0.15 0.06];
    lw = 1.25;
    fs = 8;
    fscb = 8;
    
    if nargin<3
        axi = [];
    end
    
    if isempty(axi)
        plotnew = true;
    else
        assert(length(axi)==(nr*nc))
        plotnew = false;
        ax = axi;
    end
    if nargin<4
        co = get(0,'DefaultAxesColorOrder');
        if plotnew
            color = co(1,:);
        else
            color = co(get(ax(1),'ColorOrderIndex'),:);
        end
    end


%     xl = {'$u$ (m/s)','$\psi$ (kg m/s$^2$)','T [K]','Density, $\rho_B$ (kg/m$^3$)',{'Water Mass'; 'Fraction'},{'Particle Surface'; 'Area (m^2/kg)'},{'Entrainment'; 'Coefficient, $\alpha$'}};
    xl = {'Velocity, $u$ (m/s)','Temperature, $T$ [K]','Density, $\rho$ (kg/m$^3$)',{'Water Mass'; 'Fractions'},{'Particle Surface'; 'Area, $S$ (m$^2$/kg)'},{'Entrainment'; 'Coefficient, $\alpha$'}};
   
    if plotnew
        figure('position',[50 50 1200 600])
    end
    for si = 1:(nr*nc)
        if plotnew
            ax(si) = tightSubplot(nr,nc,si,dx,dy,ppads);
        end
        
        switch si
            case 1
                p1 = plot(ax(si),wO.u,wO.z,'LineWidth',lw,'Color',color);
                if plotnew
                    ylabel(ax(si),'Height, $z$ (m a.v.l.)','Interpreter','Latex')
                end
                
%             case 2
%                 plot(ax(si),wO.psi,wO.z,'LineWidth',1.7)
                
            case 2
                plot(ax(si),wO.Tsat,wO.z,'--','LineWidth',lw,'Color',[0.4 0.4 0.4])
                hold(ax(si),'on')
%                 set(ax(si),'ColorOrderIndex',1)
                plot(ax(si),wO.T,wO.z,'LineWidth',lw,'Color',color)
                
                if plotnew
                    ty = wO.z(1) + 0.05*range(wO.z);
                    yi = find(min(abs(wO.z-ty)),1);
                    tx = wO.Tsat(yi)+15;
                    tt = text(tx,ty,'Vapor Saturation','FontSize',fs,'HorizontalAlignment','right',...
                        'VerticalAlignment','bottom','FontName','Helvetica',...
                        'Interpreter','tex','Color',[0.4 0.4 0.4]);
                    tt.Rotation = -89; %-58;
                end
                
            case 3
                plot(ax(si),wO.rho_B,wO.z,'LineWidth',lw,'Color',color)
                
            case 4
                semilogx(ax(si),wO.n_l,wO.z,'--','LineWidth',lw,'Color',color);
                hold(ax(si),'on')
                semilogx(ax(si),wO.n_v,wO.z,'LineWidth',lw,'Color',color);
%                 semilogx(ax(si),wO.n_v+wO.n_l,wO.z,'--','LineWidth',lw)
                if plotnew
                    pl(1) = plot(nan,nan,'LineWidth',lw,'Color',[0.3 0.3 0.3]);
                    pl(2) = plot(nan,nan,'--','LineWidth',lw,'Color',[0.3 0.3 0.3]);
                    [hh,icons,plots,txt] = legend(pl,{'Vapor, $n_v$','Liquid, $n_l$'},'location','southeast','FontSize',fscb,'Interpreter','Latex'); %'FontName','Helvetica');
                    icons(3).XData = [0.05 0.35];
                    icons(5).XData = [0.05 0.35];
%                     hhp = hh.Position;
%                     hh.Position = [hhp(1) hhp(2) .03 hhp(4)];
                end
                xlim(ax(si),[10^-2 1])
                set(ax(si),'XTick',10.^(-3:0))
                
            case 5
                plot(ax(si),wO.SSA,wO.z,'LineWidth',lw,'Color',color);
                
            case 6
                plot(ax(si),wO.alpha,wO.z,'LineWidth',lw,'Color',color);
        end
%         axis tight
    end
    if plotnew
        for ai = 1:length(ax)
            hold(ax(ai),'on'); 
            grid(ax(ai),'on'); 
            set(ax(ai),'FontSize',fs)
            xlabel(ax(ai),xl{ai},'Interpreter','Latex')
        end
        linkaxes(ax,'y')
        set(ax(2:end),'YTickLabel',{})
%         ylim(ax,wO.z([1 end]))
    end
    
end