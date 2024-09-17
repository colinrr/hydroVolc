function [axP,lp] = plotPlumePairMS(D,fields,z_coord,cols,linestyles,axIn)
%    [axP,lp] = plotPlumePairMS(pO,cols,fields)
%   
%   D = scalar or vector struct of hmodel outputs
%       Best to obtain D by running [D(ii).pO, D(ii).pI] = hmodel(pI);
%
% OPTIONAL IN:
%   fields      = cell vector of output fields to plot. 
%                - Nest to stack on same axes (e.g. {'T','u',{}})
%                - Will recognize and calculate: 'u_sin','n_s','n_v','n_l',
%                   and all 'atmo' props
%                - todo: custom labels for x-axes?
%   y_coord     = 'z' [Default] to use height, 'P' to use pressure levels
%   cols        = colormap matrix for custom color order
%   linestyles  = cell vector of linestyles (generally for >1 plot on same
%                   axes)
%   axIn        = replot on existing axes, must match length of fields
%
% C Rowell, Apr 2024
if nargin < 6 || isempty(axin)
    axIn = [];
end
if nargin < 5 || isempty(linestyles)
    linestyles = {'-','--'}; %get(0,'DefaultAxesLineStyleOrder');
end
if nargin < 4 || isempty(cols)
    cols = get(0, 'DefaultAxesColorOrder');
end
if nargin < 3 || isempty(z_coord)
    z_coord = 'z';
end
if nargin < 2 || isempty(fields)
    fields = {'u_sin','m_l'};
%     fields = {'rho_B','u_sin','r','theta','n_s', {'n_v','n_l'}};
end

nRuns = length(D);
nrows = 1;
ncols = length(fields);


figure('position',[50 50 1200 600])


    
    % Generate new axes as needed
    if ~isempty(axIn) && length(axIn)==ncols
        axP = axIn;
    else
        dx = 0.02;
        ppads = [0.05 0.02 0.15 0.03];

        for ai = 1:ncols
            axP(ai) = tightSubplot(nrows,ncols,ai,dx,[],ppads);
        end
        hold(axP,'on')
    end
    
    % Plot the data for each run
    for ii = 1:nRuns
        
        [zVals] = getZ(D(ii).pI,D(ii).pO,z_coord);

        ...fafo % todo: get z_tropo, z_b, z_scale properly in
            
        plotPlume(D(ii).pI, D(ii).pO, axP, fields, zVals, cols(ii,:), linestyles);
    end

    
%     if all(isfield(D,{'pI','pO'}))
%         [axP,p1] = plotPlume(D(1).pI,D(1).pO,cols(1,:));
%        
%         for ii=2:length(D)
%             p2 = plotPlume2(D(ii).pI,D(ii).pO,axP,cols(2,:));
%         end
%     end

    % Find the max height for plotting
    for ii=1:nRuns
        hm(ii) = D(ii).pO.hm;
        if z_coord == 'P'
            hmt = interpAtmoArray(D(ii).pI.atmo, hm);
            hm(ii) = hmt(1);
        end
    end
    [hm,pi] = max(hm);

    
    set(axP,'Ylim',[0 hm*zVals.ax_scale])
%     lp = [p1 p2];
    linkaxes(axP,'y')
    
    

end

function [zVals] = getZ(pI, pO, z_coord)
% Return either z (m) or P (Pa) for vertical axis
%  + z_tropo (tropopause height), z_b (neutral buoyancy), ax_scale (unit
%  conversion from standard SI)

h_scale = 1/1e3; % m to km
p_scale = 1/100; % Pa to hPa

switch z_coord
    case 'z'
        z = pO.z .* h_scale;
        
        if ~isempty(pO.atmo.ztropo)
            z_tropo = pO.atmo.ztropo .* h_scale;
        else
            z_tropo = nan;
        end
        z_b = pO.hb .* h_scale;
        
        ax_scale = h_scale;
        z_unit = 'km';
        
    case 'P'
        interpMethod = 'pchip';
        z = 10.^interp1(pI.atmo(:,2),log10(pI.atmo(:,1)),real(pO.z),interpMethod,'extrap') .* p_scale;
        
        if ~isempty(pO.atmo.ztropo)
            z_tropo = 10.^interp1(pI.atmo(:,2),log10(pI.atmo(:,1)),real(pO.atmo.ztropo),interpMethod,'extrap') .* p_scale;
        else
            z_tropo = nan;
        end        
        z_b = 10.^interp1(pI.atmo(:,2),log10(pI.atmo(:,1)),real(pO.hb),interpMethod,'extrap') .* p_scale;
        ax_scale = p_scale;
        z_unit = 'hPa';
        
    otherwise
        error('Z coordinate choice not recognized')
end

zVals.z = z;
zVals.z_tropo = z_tropo;
zVals.z_b = z_b;
zVals.ax_scale = ax_scale;
zVals.unit = z_unit;

end

function [axP,lp] = plotPlume(pI, pO, axP, fields, z, color, linestyles)

    lw = 1.25; % Line Width  
    fs = 8; % Font size
    fscb = 8; %Legend font size
    
    ptp = ~isempty(pO.atmo.ztropo);
    
    ncol = length(fields);
    
    defaultPlots(pI, pO, axP, z, color)
    
%  % FILL IN PROPER BEHAVIOR LATER (or just get into python)     
%     for fi = 1:ncol
%         
%         % Field recognition here
%         if iscell(fields{fi})
%             for ji=1:length(fields{fi})
%                 plot(axP(fi), z, pO.(fields{fi}),'LineWidth',lw,'color',color,'linestyle',linestyles{ji})
%             end
%             
%         else
%             plot(axP(fi), z, pO.(fields{fi}),'LineWidth',lw,'color',color)
%         end
%         
%         % Add tropopause
%         if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw); end
%         
%         % Add LNB
%         
%     end
end

function [plotVars, axlabels, leglabels] = getPlotVars(pO,fields)
% Lookup function here for conversions, labels, units?
 %  atmo fields get assigned black color?
plotVars = cell(1,length(fields));

for fi=1:length(fields)
    if ischar(fields{fi})
        fname = fields(fi);
        
        switch fname
            case 'u_sin' %'u_sin','n_s','n_v','n_l','atmo' props
                plotVars{fi} = pO.u.*sin(pO.angle);
                labels{fi} = 'Vert. Velocity (m/s)';
                
            case 'n_s'
                plotVars{fi} = (pO.m_s)./pO.m;
                labels{fi} = {'Particle Mass'; 'Fraction, $n_s$'};

            case 'n_v'
                plotVars{fi} = (pO.m_v)./pO.m;
                labels{fi} = 'Vapor, $n_v$';
                
            case 'n_l'
                plotVars{fi} = (pO.m_l)./pO.m;
                labels{fi} = 'Liquid, $n_l$';
                
%             case 'theta_a'
        end
        
    elseif iscell(fields{fi})
        [plotVars{fi},labels{fi}] = getPlotVars(pO,fields{fi});
    end
    
    
    
end
end

function defaultPlots(pI, pO, axP, z, color)
% Temporary easy thing to just fill in needed fields

    lw = 1.25; % Line Width  
    fs = 8; % Font size
    fscb = 8; %Legend font size


% Velocity
    plot(axP(1), pO.u.*sin(pO.angle),z.z,...
    'LineWidth',lw,'Color',color(1,:))
    
    % Ambient
%     ax(np) = gca;
    ylim([0 max(pO.z)])
    xl = [0 200]; %xlim;
    plot(axP(1),xl,z.z_tropo*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw)
    plot(axP(1),xl,z.z_b*[1 1],':','LineWidth',lw,'Color',color(1,:))
    axP(1).TickLabelInterpreter = 'tex';
    xlabel({'Velocity, $u$ (m/s)'},'Interpreter','Latex')
    axP(1).FontSize = fs;
    grid(axP(1),'on')    
    
% Saturation
%     epsilon = 287/461;
% %     pVals = getZ(pI, pO, 'P');
%     P_a = 10.^interp1(pI.atmo(:,2),log10(pI.atmo(:,1)),real(pO.z),'pchip','extrap');
%     es   = 100.*6.112.*exp(17.67.*(pO.theta - 273.15)./(pO.theta +243.5-273.15)); % saturation vapor pressure for water
%     dp = max([P_a-es  2*es],[],2); % Avoid a div0 when abs. pressure becomes very small
% %     w_s1    = 1./epsilon .* es ./ dp; % w_a = r_h*w_s, mass mixing ratio of water vapor to dry air at saturation
%     w_s    = 1./epsilon .* es ./ (P_a-es);
%     w_s(es>=P_a)=nan; % Nan out the pure vapour regime for now
%     n_v    = pO.m_v ./ pO.m_d; % Actual vapour mass mixing ratio;
%     
%     plot(axP(2), n_v,z.z,...
%     'LineWidth',lw,'Color',color(1,:))
%     xl = [0 2];
%     axP(2).TickLabelInterpreter = 'tex';
%     xlabel(axP(2),{'Vapor'; 'Saturation'},'Interpreter','Latex')
%     axP(2).FontSize = fs;
%     grid(axP(2),'on')    
%     plot(axP(2),xl,z.z_tropo*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw)
%     plot(axP(2),xl,z.z_b*[1 1],':','LineWidth',lw,'Color',color(1,:))
    
% Liquid mass
    plot(axP(2),pO.m_l/(pO.m_l + pO.m_v),z.z,...  % liquid mass fraction
    'LineWidth',lw,'Color',color)
    set(axP(2),'XScale','log')
    xl = xlim(axP(2));
    plot(axP(2),xl,z.z_tropo*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',lw)
    plot(axP(2),xl,z.z_b*[1 1],':','LineWidth',lw,'Color',color(1,:))
    axP(2).FontSize = fs;
    grid(axP(2),'on')    
    xlabel(axP(2),{'Liquid water'; 'mass (kg)'},'Interpreter','Latex')


end