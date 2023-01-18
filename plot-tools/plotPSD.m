function plotPSD(pI,usePhi)
 % Show grain size distribution parameters
 % pI = plume input struct (see getPlumeSource)
 % usePhi = plot with phi size units. False uses particle radii
 
 if nargin<2
     usePhi = true;
 end
 
 % Mass fraction
 % Cumulative # (power law equiv)
 % Porosity or gas mass fraction
 % Specific surface area?
 
    nr = 4;
    nc = 1;
    dx = [];
    dy = 0.01;
    ppads = [0.1 0.03 0.12 0.05];
    lw = 1.7;
    fs = 12;
 
    phi = -log2(2000*pI.Rgsd);
    Ni = 2.^(log2(1) + pI.D.*phi);
    Ni_all = cumsum(Ni./Ni(1));  % Cumulative N (relative)
    SSA = (3*pI.saScale*sum(pI.Rgsd.^(2-pI.D))./sum(pI.rhoi.*pI.Rgsd.^(3-pI.D))); % Total SSA
    SSAi = 3*pI.saScale./(pI.rhoi.*pI.Rgsd).*pI.nsi; % Specific surface area per particle bin per bulk mass
%     SSAi = SSA.*pI.nsi; % Specific surface area per particle bin
    
    if usePhi
        x = phi;
        xl = '\phi';
        useLogx = false;
    else
        x = pI.Rgsd*2/1e3;
        xl = 'd (mm)';
        useLogx = true;
    end
    yl = {'Mass %','Cumulative N','Vesicle Mass %','Particle SSA (m^2/kg)'};
    
    figure('position',[50 50 600 750])   
    for si = 1:(nr*nc)
        ax(si) = tightSubplot(nr,nc,si,dx,dy,ppads);
        co = get(gca,'ColorOrder');
        co1 = co(1,:);
        
        switch si
            case 1
%                 plot(ax(si),x,pI.nsi,'LineWidth',1.7)
                b1=bar(ax(si),x,pI.nsi,'histc');
                set(b1,'FaceColor','None','EdgeColor',co1,'LineWidth',1.5)
                ylim(ax(si),[0 max(pI.nsi)*1.1])
                set(gca,'XTickLabel',[])
                
            case 2
                semilogy(ax(si),x,Ni_all,'LineWidth',1.7)
                set(gca,'XTickLabel',[])

            case 3
                b2=bar(ax(si),x,pI.ni,'histc');
                set(b2,'FaceColor','None','EdgeColor',co1,'LineWidth',1.5)
                ylim(ax(si),[0 max(pI.ni)*1.1])
                set(gca,'XTickLabel',[])
                
            case 4
                b3=bar(ax(si),x,SSAi,'histc');
                set(b3,'FaceColor','None','EdgeColor',co1,'LineWidth',1.5)
                ll = sprintf('SA Scale = %.0f\nSSA = %.3e m^2/kg',pI.saScale,SSA);
                text(ax(si),0.05, 0.8,ll,'Units','Normalized','FontSize',fs,'HorizontalAlignment','left')
                ylim(ax(si),[0 max(SSAi)*1.1])
                xlabel(xl)
                
        end
        ylabel(yl{si})
%         axis tight

    end
    set(ax,'FontSize',fs)
    hold(ax,'on')
    linkaxes(ax,'x')
    
end