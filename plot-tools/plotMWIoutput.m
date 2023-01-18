function [ax,p1] = plotMWIoutput(D,psd,diag,cols)
% plotMWIouput(D)
% plots the detailed output of a single coupled conduit-plume run
% D = struct of coupled model output. For struct array, will plot each in
%   sequence on the same figure
% psd = flag to plot psd
% diag = flag to plot integration stepping
%
% CRowell Jun 2022

    if nargin<4
        cols = get(0,'DefaultAxesColorOrder');
%         cols = cols([2 1 3 4 5 6 7],:);
    end
    if nargin<3
        diag = false;
    end
    if nargin<2
        psd = true;
    end
    
%     nrows = 1;
%     ncols = 6;
    
    for ii=1:length(D)
        zm(ii) = D(ii).cI.Zw - D(ii).wO.z(1);
    end
    [zm,pi] = max(zm);


    figure('position',[50 50 1350 600])
    
    if all(isfield(D,{'pI','pO'}))
        [axP,p1] = plotMWI(D(1).wO,D(1).cI,D(1).pI,[],cols(1,:),psd,diag);
       
        for ii=1:length(D)
%             if ii~=pi
                p2 = plotMWI(D(ii).wO,D(ii).cI,D(ii).pI,axP,cols(ii,:),psd,diag);
%             end
        end
    end
    set(axP,'Ylim',[-zm 0])
    lp = [p1 p2];

end


function [ax,p1] = plotMWI(wO,cI,pI,axi,color,psd,diag)

    % U, psi, T, rho, n_l, n_v, A
    nr = 1;

    dx = 0.02;
    dy = [];
    ppads = [0.06 0.02 0.15 0.06];
    lw = 1.25;
    fs = 8;
    fscb = 8;
    
    PSDdz = 5;
    
    if nargin<7
        diag = false;
    end
%     if nargin<6
    if xor(diag,psd)
        nc = 7;
        if diag
            psdAx = [];
            diagAx = 7;
        else
            psdAx = 7;
            diagAx = [];
        end
    elseif and(diag,psd)
        nc = 8;
        psdAx = 7;
        diagAx = 8;
    elseif and(~diag,~psd)
        nc = 6;
        psdAx = [];
        diagAx = [];
    end
    
    if nargin<4
        axi = [];
    end    
    if isempty(axi)
        plotnew = true;
    else
        assert(length(axi)==(nr*nc))
        plotnew = false;
        ax = axi;
    end
    if nargin<5
        co = get(0,'DefaultAxesColorOrder');
        if plotnew
            color = co(1,:);
        else
            color = co(get(ax(1),'ColorOrderIndex'),:);
        end
    end

%     xl = {'$u$ (m/s)','$\psi$ (kg m/s$^2$)','T [K]','Density, $\rho_B$ (kg/m$^3$)',{'Water Mass'; 'Fraction'},{'Particle Surface'; 'Area (m^2/kg)'},{'Entrainment'; 'Coefficient, $\alpha$'}};
    xl = {'Velocity, $u$ (m/s)',{'Entrainment'; 'Coefficient, $\alpha$'},...
        'Temperature, $T$ [K]','Density, $\rho$ (kg/m$^3$)',{'Water Mass'; 'Fractions'},...
        {'Particle Surface'; 'Area, $S$ (m$^2$/kg)'},{'PSD ($\phi$)'},'steps/m'};
    if ~diag
        xl(8) = [];
    end
    if ~psd
        xl(7) = [];
    end
   
%     if plotnew
%         figure('position',[50 50 1200 600])
%     end
    for si = 1:(nr*nc)
        if plotnew
            ax(si) = tightSubplot(nr,nc,si,dx,dy,ppads);
        end
        z = wO.z - cI.Zw;
        
        switch si
            case 1
                p1 = plot(ax(si),wO.u,z,'LineWidth',lw,'Color',color);
                if plotnew
                    ylabel(ax(si),'Depth, $z$ (m b.w.l.)','Interpreter','Latex')
                end
                
%             case 2
%                 plot(ax(si),wO.psi,wO.z,'LineWidth',1.7)
                
            case 2
                plot(ax(si),wO.alpha,z,'LineWidth',lw,'Color',color);

            case 3
                if plotnew
                    plot(ax(si),wO.Tsat,z,'--','LineWidth',lw,'Color',[0.4 0.4 0.4])
                    hold(ax(si),'on')
                end
%                 set(ax(si),'ColorOrderIndex',1)

                tgi1 = find(wO.T<=pI.T_g+pI.T_g_rng,1,'first');
                tgi2 = find(wO.T>=pI.T_g,1,'last');
                if and(~isempty(tgi1),~isempty(tgi2))
                    tgi = sort([tgi1 tgi2]);
%                     tgi = [tgi fliplr(tgi)];
                    cc = rgba2rgb(color,0.4);
                    patch(ax(si),wO.T(tgi([1 2 2 1])),z(tgi([1 1 2 2])),cc,'EdgeAlpha',0 )
                end
                plot(ax(si),pI.T_g.*[1 1],z([1 end]),'--','LineWidth',lw,'Color',color)
                
                plot(ax(si),wO.T,z,'LineWidth',lw,'Color',color)
                
                if plotnew
                    ty = wO.z(1) + 0.05*range(wO.z);
                    yi = find(min(abs(wO.z-ty)),1);
                    tx = wO.Tsat(yi)+15;
                    tt = text(tx,ty,'Vapor Saturation','FontSize',fs,'HorizontalAlignment','right',...
                        'VerticalAlignment','bottom','FontName','Helvetica',...
                        'Interpreter','tex','Color',[0.4 0.4 0.4]);
                    tt.Rotation = -89; %-58;
                end
                
            case 4
                plot(ax(si),wO.rho_B,z,'LineWidth',lw,'Color',color)
                
            case 5
                semilogx(ax(si),wO.n_l,z,'--','LineWidth',lw,'Color',color);
                hold(ax(si),'on')
                semilogx(ax(si),wO.n_v,z,'LineWidth',lw,'Color',color);
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
                
            case 6
                if and(~isempty(tgi1),~isempty(tgi2))
                    tgi = sort([tgi1 tgi2]);
%                     tgi = [tgi fliplr(tgi)];
                    cc = rgba2rgb(color,0.4);
                    patch(ax(si),wO.SSA(tgi([1 2 2 1])),z(tgi([1 1 2 2])),cc,'EdgeAlpha',0 )
                    hold(ax(si),'on')
                end
                plot(ax(si),wO.SSA,z,'LineWidth',lw,'Color',color);

            case psdAx
                PSDz = sort(0:-PSDdz:-cI.Zw)';
                zi = and( PSDz>=wO.z(1)-cI.Zw, PSDz<=(wO.z(max(tgi))-cI.Zw) );
                [nsi_scaled,~] = makePSDplot(wO,PSDz(zi)+cI.Zw);
                plot(ax(si),wO.phi,nsi_scaled-cI.Zw,'LineWidth',lw,'Color',color);
                if plotnew
                    set(ax(si),'YTick',PSDz)
                end
                
            case diagAx
                zdens = [0; 1./diff(wO.z)];
                semilogx(ax(si),zdens,z,'LineWidth',lw,'Color',color);
                
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

function [nsi,zz] = makePSDplot(wO,zz)
    zscale = 10;
%     zz = linspace(min(zlim),max(zlim),nn)';
    dz = mean(diff(zz));
    zi = zeros(length(zz),1);
    for ii = 1:length(zz)
        [~,zi(ii)] = min(abs(wO.z-zz(ii)));
    end
    nsi = (wO.nsi(zi,:).*dz.*zscale + zz)';
%     zz  = zz(zi);
end

function rgb = rgba2rgb(rgb,Alpha,rgb_bg)
% rgb = rgba2rgb(rgb,Alpha,rgb_bg)
% Convert rgb with transparency to an equivalent rgb.
%  rgb = original rgb value, range (0 1)
%  alpha = transparency value (0 1)
%  rgb_bg = background rgb (defaults to white, [1 1 1])
%
%   rgb can be a single 1x3 color vector or nx3 color matrix
%   alpha can be a single scaler or vector of length == n
%
%     return new Color(
%         (1 - alpha) * RGB_background.r + alpha * RGBA_color.r,
%         (1 - alpha) * RGB_background.g + alpha * RGBA_color.g,
%         (1 - alpha) * RGB_background.b + alpha * RGBA_color.b

    if nargin<3
        rgb_bg = [1 1 1];
    end
    
    numcolors = size(rgb,1);
    numalph = numel(Alpha);
    nn = max([numcolors numalph]);
    assert(or(numalph==1,numalph==nn),'Alpha must be a scalar or vector with length = size(rgb,1)')

    for ii=1:nn
        if numalph==nn
            alpha = Alpha(ii);
        elseif numalph==1
            alpha = Alpha;
        end

        if numcolors ==nn
            r = rgb(ii,1);
            g = rgb(ii,2);
            b = rgb(ii,3);
        elseif numcolors==1
            r = rgb(1);
            g = rgb(2);
            b = rgb(3);
        end

        r2 = (1 - alpha) * rgb_bg(1) + alpha * r;
        g2 = (1 - alpha) * rgb_bg(2) + alpha * g;
        b2 = (1 - alpha) * rgb_bg(3) + alpha * b;

        rgb(ii,:) = [r2 g2 b2];
    end
end