function plotCoupledPlume(d)
% plotCouplePlume(D)
% plots the detailed output of a single coupled conduit-plume run
%
% CRowell Mar 2021

    
    nrows = 2;
    ncols = 6;

    figure('position',[50 50 1200 1200])

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

    
    if all(isfield(d,{'conduitI','conduitO'}))
        Conduit_flow_plot(d.conduitI,d.conduitO)
        
        if d.conduitO.Par.frag
            plotPlume(d.plumeI,d.plumeO)
        end
    elseif all(isfield(d,{'plumeI','plumeO'}))
        plotPlume(d.plumeI,d.plumeO)
       
    end
    
end


function Conduit_flow_plot(varargin)
% Input either single path to results file, or two args: Input, Output 
% structs from conduit flow model
% Borrowed from Hajimirza code

% clear 
% close all

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);
set(0, 'DefaultTextFontSize', 7);
set(0, 'DefaultPatchEdgeColor','none');
set(0, 'DefaultTextInterpreter', 'tex');
set(0,'DefaultAxesBox','on');


if nargin==1 && ischar(varargin{1})
    load(varargin{1})
elseif nargin==1
    error('Single input must be a valid .mat results file')
elseif nargin==2
    Output = varargin{2};
    Input  = varargin{1};
    assert(and(istable(Input),isstruct(Output)),'Input must be table, Output must be struct')
end



frag_index = find(Output.porosity>=Input.phi_frag,1);


% ylim
yl = 6;
color = [0 0 .8];
    
lw = 1; % Line Width  
fs = 10; % Font size
    

% figure('position',[50 50 1200 1200])
nf = 6; % Number of subplots
nc = 2;
np = 6;

%===================================================
% Pressure
np = np + 1;
subplot(nc,nf,np)
hold on

   
plot(Output.pm/1e6,...
     Output.Z/1e3,...
    'LineWidth',lw,'Color',color(1,:));        


% Lithostatic pressure
Z = linspace(0,10000,1000);
plot(Z*2400*9.81/1e6,Z/1e3,'k--');
  


ax(np) = gca;
ax(np).YDir = 'reverse';
xlim([0 150])
ylim([0 yl])
ax(np).TickLabelInterpreter = 'tex';
xlabel({'Pressure (MPa)'})
ax(np).XTick = [0:50:400];
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;
ylabel('Depth, Z (km)')
grid on

text(ax(np).XLim(2)-10,.3,...
    '\bf{(a)}','FontSize',fs,'FontWeight','bold',...
    'VerticalAlignment', 'bottom') 

tt = text(500*2400*9.81/1e6+10,.5,'Lithostatic','FontSize',fs);
tt.Rotation = -58;


%===================================================
% U
np = np + 1;
subplot(nc,nf,np)
hold on



pp1 = plot(Output.U,...
           Output.Z/1e3,...
           'LineWidth',lw,'Color',color(1,:));   
    
pp2 = plot(Output.U./Output.M,...
           Output.Z/1e3,':',...
           'LineWidth',lw,'Color',color(1,:));    




ax(np) = gca;
ax(np).YDir = 'reverse';
xlim([0 1000])
ylim([0 yl])
ax(np).TickLabelInterpreter = 'tex';
xlabel({'Velocity (m/s)'})
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;

grid on


text(ax(np).XLim(1)+10,.3,...
    '\bf{(b)}','FontSize',fs,'FontWeight','bold',...
    'VerticalAlignment', 'bottom') 


legend([pp1,pp2],'magma','sound','Location','southeast')

%===================================================
% Conduit
np = np + 1;
subplot(nc,nf,np)
hold on

   
plot(Output.dPdt/1e6,...
     Output.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
ax(np).XScale = 'log';
xlim([1e-2 1e2])
ylim([0 yl])
ax(np).XTick = 10.^([-4:2:2]);
ax(np).TickLabelInterpreter = 'tex';
xlabel({'Decompression rate' '(MPa s^{-1})'})
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;



text(ax(np).XLim(2)-40,...
    .3,...
    '\bf{(c)}','FontSize',fs,'FontWeight','bold',...
    'VerticalAlignment', 'bottom') 

% grid on



%===================================================
% Gas volume fraction
np = np + 1;
subplot(nc,nf,np)
hold on
   
plot(Output.porosity * 100,...
     Output.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
xlim([0 100])
ylim([0 yl])
ax(np).XTick = [0:20:100];
ax(np).TickLabelInterpreter = 'tex';
xlabel({'Gas volume fraction (%)'})
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;
ylabel('Depth, Z (km)')

grid on

text(ax(np).XLim(1),.3,...
    '\bf{(d)}','FontSize',fs,'FontWeight','bold',...
    'VerticalAlignment', 'bottom') 



%===================================================
% Nucleation rate
np = np + 1;
subplot(nc,nf,np)
hold on



Output.J(Output.J<1e-10) = 1e-10;
plot(Output.J,...
     Output.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
ax(np).XScale = 'log';
% xlim([1e6 1e16])
ylim([0 yl])
% ax(np).XTick = 10.^([6:4:20]);
ax(np).TickLabelInterpreter = 'tex';
xlabel({'Nucleation rate (m^{-3} s^{-1})'})
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;


% grid on

plot([1e-3 1e-3],[2.65,.8],'Color','k');



 text(ax(np).XLim(1),.3,...
    '\bf{(e)}','FontSize',fs,'FontWeight','bold',...
    'VerticalAlignment', 'bottom') 


%===================================================
% BND
np = np + 1;
subplot(nc,nf,np)
hold on


plot(Output.M0,...
     Output.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
ax(np).XScale = 'log';
% xlim([1e8 1e16])
ylim([0 yl])
% ax(np).XTick = 10.^([2:2:31]);
ax(np).TickLabelInterpreter = 'tex';
xlabel({'Bubble number density (m^{-3})'})
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;


text(ax(np).XLim(1),.3,...
    '\bf{(f)}','FontSize',fs,'FontWeight','bold',...
    'VerticalAlignment', 'bottom') 




%========================================================
% Setting up subplot positions
dx1 = .05;
x2 = .27;
dx2 = .03;
x1 = x2;

np = 1+nf;
% ax(np).Position = [dx1 .15 x1 .8];


% for np = (1+nf):(nf*nc)
%     ax(np).Position = [dx1+x1+(np-2)*(x2+dx2)+dx2 .55 x2 .40];
%     ax(np+nf).Position = [dx1+x1+(np-2)*(x2+dx2)+dx2 .07 x2 .40];
% end

for np = (2+nf):(nf*nc)
    ax(np).YTickLabel = {};
    ax(np).YTick = [0:20];
%     ax(np+nf).YTickLabel = {};
%     ax(np+nf).YTick = [0:20];

end

%========================================================
set(gcf,'PaperUnits','centimeters')
figuresize = [18 18];
set(gcf,'PaperPosition', [0 0.1 figuresize(1) figuresize(2)]);
set(gcf,'PaperSize', figuresize);
% print('Results','-dpdf')




end

function plotPlume(pI,pO)



    % ylim
    yl = 6;
    color = [0 0 .8];

    lw = 1; % Line Width  
    fs = 10; % Font size


%     figure
    nf = 6; % Number of subplots
    nc = 2;
    np = 0;
    
    ptp = ~isempty(pO.atmo.ztropo);
    
    % ===================================================
    % Density
    np = np + 1;
    subplot(nc,nf,np)
    hold on
    
    plot(pO.rho_B,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
    
    % Ambient
    yl = [0 max(pO.z)/1e3];
    plot(pO.atmo.rho,pO.z/1e3,'--k')
    xl = xlim;
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(xl,pO.hb/1e3*[1 1],'--r')
    ax(np) = gca;
%     xlim([0 150])
    ylim(yl)
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Density [kg/m^3]'})
%     ax(np).XTick = [0:50:400];
    ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    ylabel('Height, Z [km a.v.l.]')
    grid on
    
    % ===================================================
    % Velocity
    np = np + 1;
    subplot(nc,nf,np)
    hold on
    
    plot(pO.u,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
    
    % Ambient
%     plot(pO.rho_aB,pO.z/1e3,'--k')
    ax(np) = gca;
%     xlim([0 150])
    ylim([0 max(pO.z)/1e3])
    xl = xlim;
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(xl,pO.hb/1e3*[1 1],'--r')
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Velocity, U [m/s]'})
%     ax(np).XTick = [0:50:400];
    ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on

    % ===================================================
    % Radius
    np = np + 1;
    subplot(nc,nf,np)
    hold on
    
    plot(pO.r,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
%     set(gca,'XScale','log') % YEH??

    % Ambient
%     plot(pO.rho_aB,pO.z/1e3,'--k')
    xl = [1 1e4];
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(xl,pO.hb/1e3*[1 1],'--r')

    ax(np) = gca;
    xlim(xl)
    ylim([0 max(pO.z)/1e3])
%     ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Radius [m]'})
%     ax(np).XTick = [0:50:400];
    ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on    

    % ===================================================
    % Temperature
    np = np + 1;
    subplot(nc,nf,np)
    hold on
    
    plot(pO.theta,pO.z/1e3,...
    'LineWidth',lw,'Color',color(1,:))
    
    % Ambient
    plot(pO.atmo.theta_a,pO.z/1e3,'--k')
    ax(np) = gca;
    xl = xlim;;
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(xl,pO.hb/1e3*[1 1],'--r')

    %     xlim([0 150])
    ylim([0 max(pO.z)/1e3])
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Temperature [K]'})
%     ax(np).XTick = [0:50:400];
    ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on   
    
    % ===================================================
    % Mass or volume fractions?
    np = np + 1;
    subplot(nc,nf,np)
    hold on
    
    % Ambient vapour mixing ratio
    plot(pO.atmo.w_a,pO.z/1e3,'--',...
    'LineWidth',lw,'Color',color(1,:))
    % Saturation vapour mixing ratio 
    plot(pO.atmo.w_s,pO.z/1e3,'--k',...
    'LineWidth',lw)

    plot((pO.m-pO.m_d-pO.m_l-pO.m_v)./pO.m,pO.z/1e3,...  % particle mass fraction
    'LineWidth',lw)
    plot((pO.m_v)./pO.m,pO.z/1e3,...  % vapour mass fraction
    'LineWidth',lw)
    plot((pO.m_d)./pO.m,pO.z/1e3,...  % dry air mass fraction
    'LineWidth',lw)    
    plot((pO.m_l)./pO.m,pO.z/1e3,...  % liguid water mass fraction
    'LineWidth',lw)    
    set(gca,'XScale','log')
    
    xl = [1e-7 1e0];
    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(xl,pO.hb/1e3*[1 1],'--r')

    ax(np) = gca;
    xlim(xl)
    ylim([0 max(pO.z)/1e3])
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Mass fractions'})
%     ax(np).XTick = [0:50:400];
    ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on   
    
    % ===================================================
    % GSD
    np = np + 1;
    subplot(nc,nf,np)
    hold on
    

    pcolor(log10(pI.Rgsd),pO.z/1e3,pO.m_si)
    shading flat
    colormap(gray)
    xl = xlim;

    if ptp; plot(xl,pO.atmo.ztropo/1e3*[1 1],'--','Color',[0.5 0.5 0.5]); end
    plot(xl,pO.hb/1e3*[1 1],'--r')

    ax(np) = gca;
    xlim(xl)
    ylim([0 max(pO.z)/1e3])
    ax(np).TickLabelInterpreter = 'tex';
    xlabel({'Mass fractions'})
%     ax(np).XTick = [0:50:400];
    ax(np).FontName = 'Arial';
    ax(np).FontSize = fs;
    grid on   

end