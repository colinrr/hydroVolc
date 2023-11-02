function varargout = plotConduitOutput(D,cols)
% plotCouplePlume(D)
% plots the detailed output of a pair of conduit runs
% d1 = struct containing conduit model input/output (may be struct vector):
%        D.cI = output of getConduitSource
%        D.cO = output of Conduit_flow_with_nucleation

% Same for d2
% cols = 2x3 color vector for the 2 plots
%
% CRowell Aug 2021

    if nargin<2
        if length(D)<=6
            cols = get(0,'DefaultAxesColorOrder');
            cols = cols([1:2 4:end],:);
        else
            cols = viridis(length(D));
        end 
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

    zm = zeros(length(D),1);
    for ii=1:length(D)
        zm(ii) = D(ii).cI.Z0;
    end
    [zm,~] = max(zm);

    
    if all(isfield(D,{'cI','cO'}))
        
        
        if length(D)>1
            [axC,lh(length(D))] = Conduit_flow_plot(D(end).cI,D(end).cO,cols(end,:));
            for ii=length(D):-1:1
                lh(ii) = Conduit_plot_2(D(ii).cI,D(ii).cO,axC,cols(ii,:));
                
            end
        else
            [axC,lh] = Conduit_flow_plot(D.cI,D.cO,cols(1,:));
        end
    end
    
    ylim(axC,[0 zm/1e3])
    linkaxes(axC,'y')
    
    if nargout>=1
        varargout{1} = axC;
    end
    if nargout==2
        varargout{2} = lh;
    end
    

end

function [ax,p1] = Conduit_flow_plot(cI,cO,color)
% Input either single path to results file, or two args: Input, Output 
% structs from conduit flow model
% Borrowed from Hajimirza code

% clear 
% close all

% set(0, 'DefaultAxesFontName', 'Arial');
% set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 9);
set(0, 'DefaultTextFontSize', 9);
set(0, 'DefaultPatchEdgeColor','none');
set(0, 'DefaultTextInterpreter', 'Latex');
set(0,'DefaultAxesBox','on');

dx = 0.02;
ppads = [0.05 0.02 0.15 0.05];

% if nargin==1 && ischar(varargin{1})
%     load(varargin{1})
% elseif nargin==1
%     error('Single input must be a valid .mat results file')
% elseif nargin==2
%     Output = varargin{2};
%     Input  = varargin{1};
% %     assert(and(istable(Input),isstruct(Output)),'Input must be table, Output must be struct')
% end



frag_index = find(cO.porosity>=cI.phi_frag,1);


% ylim
yl = 6;
% color = [0 0 .8];
    
lw = 1.25; % Line Width 
lw2 = 1;
fs = 8; % Font size
fscb = 8;

% figure('position',[50 50 1200 1200])
nf = 6; % Number of subplots
nc = 1;
np = 0;

%===================================================
% Pressure
np = np + 1;
ax(1) = tightSubplot(nc,nf,np,dx,[],ppads);
hold on

% Lithostatic pressure
plith = (cO.Z*2400*9.81+cI.pf);
  
% Straight up pressures
p1 = plot(cO.pm/1e6,...
     cO.Z/1e3,...
    'LineWidth',lw,'Color',color(1,:));        
plot(plith/1e6,cO.Z/1e3,'--','Color',color(1,:));
xlab = {'Pressure (MPa)'};  


% Overpressure ratio
% plot(cO.pm./plith,...
%      cO.Z/1e3,...
%     'LineWidth',lw,'Color',color);        
% xlab = '$\frac{p}{p_{lith}}$';


ax(np) = gca;
ax(np).YDir = 'reverse';
% xlim([0 150])
ylim([0 yl])
% ax(np).TickLabelInterpreter = 'tex';
xlabel(xlab)
% ax(np).XTick = [0:50:400];
% ax(np).FontName = 'Arial';
ax(np).FontSize = fs;
ylabel('Depth, $z$ (km)')
grid on
hold on
xl = xlim;
plot(xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  

% text(ax(np).XLim(2)-10,.3,...
%     '\bf{(a)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 

tt = text(3100*2400*9.81/1e6+2,3.1,'Lithostatic','FontSize',fs,'FontName','Helvetica',...
    'Interpreter','tex','VerticalAlignment','bottom');
tt.Rotation = -72; %-58;


%===================================================
% U or M
np = np + 1;
ax(2) = tightSubplot(nc,nf,np,dx,[],ppads);
hold on

% pp1 = plot(cO.U,...
%            cO.Z/1e3,...
%            'LineWidth',lw,'Color',color(1,:));   
    
pp2 = plot(cO.M,...
           cO.Z/1e3,...
           'LineWidth',lw,'Color',color(1,:));    

ax(np) = gca;
ax(np).YDir = 'reverse';
% xlim([0 1000])
xlim([0 1])
ylim([0 yl])
% ax(np).TickLabelInterpreter = 'tex';
% xlabel({'Velocity (m/s)'})
xlabel('Mach Number')
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;
hold on
xl = xlim;
plot(xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  
tf = text(0.98,cO.Par.Zf/1e3-0.15,'Fragmentation','FontSize',fs,'FontName','Helvetica',...
    'Interpreter','tex','VerticalAlignment','bottom','HorizontalAlignment','right');

grid on

% text(ax(np).XLim(1)+10,.3,...
%     '\bf{(b)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 


% legend([pp1,pp2],'magma','sound','Location','southeast')

%===================================================
% Decompression rate
np = np + 1;
ax(3) = tightSubplot(nc,nf,np,dx,[],ppads);
hold on

   
plot(cO.dPdt/1e6,...
     cO.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
ax(np).XScale = 'log';
xlim([1e-2 1e1])
ylim([0 yl])
ax(np).XTick = 10.^([-2:1]);
ax(np).XTickLabel = {'10^{-2}','10^{-1}','10^0','10^{1}'}; %10.^([-2:16]);
% ax(np).XTickLabel = {'10^{-2}','','10^0','','10^2'}; %10.^([-2:16]);
% ax(np).TickLabelInterpreter = 'tex';
xlabel({'Decompression'; 'Rate (MPa s$^{-1}$)'})
% ax(np).FontName = 'Arial';
ax(np).FontSize = fs;
hold on
xl = xlim;
plot(xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  

grid on 
set(ax(np),'XMinorTick','off')
set(ax(np),'YMinorTick','off')
% text(ax(np).XLim(2)-40,...
%     .3,...
%     '\bf{(c)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 

% grid on



%===================================================
% Gas volume fraction
np = np + 1;
ax(4) = tightSubplot(nc,nf,np,dx,[],ppads);
hold on
   
plot(cO.porosity * 100,...
     cO.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
xlim([0 100])
ylim([0 yl])
% ax(np).XTick = [0:20:100];
% ax(np).TickLabelInterpreter = 'tex';
xlabel({'Gas Volume'; 'Fraction (\%)'})
% ax(np).FontName = 'Arial';
ax(np).FontSize = fs;
% ylabel('Depth, Z (km)')
hold on
xl = xlim;
plot(xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  

grid on

% text(ax(np).XLim(1),.3,...
%     '\bf{(d)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 



%===================================================
% Nucleation rate
np = np + 1;
ax(5) = tightSubplot(nc,nf,np,dx,[],ppads);
hold on



cO.J(cO.J<1e-10) = 1e-10;
plot(cO.J,...
     cO.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
ax(np).XScale = 'log';
xlim([1e6 1e16])
ylim([0 yl])
ax(np).XTick = 10.^([6:4:20]);
% ax(np).TickLabelInterpreter = 'tex';
xlabel({'Nucleation'; 'Rate (m$^{-3}$ s$^{-1}$)'})
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;
hold on
xl = xlim;
plot(xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  


grid on 
set(ax(np),'XMinorTick','off')
set(ax(np),'YMinorTick','off')


plot([1e-3 1e-3],[2.65,.8],'Color','k');



%  text(ax(np).XLim(1),.3,...
%     '\bf{(e)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 


%===================================================
% Supersaturation
np = np + 1;
ax(6) = tightSubplot(nc,nf,np,dx,[],ppads);
hold on


% plot(cO.M0,...
%      cO.Z/1e3,...
%      'LineWidth',lw,'Color',color(1,:));        
% plot((cO.Cm-cO.Csol)*100,cO.Z/1e3,'LineWidth',lw,'Color',color(1,:)); 
plot((cO.Cm)*100,cO.Z/1e3,'LineWidth',lw,'Color',color(1,:)); 

ax(np) = gca;
ax(np).YDir = 'reverse';
% ax(np).XScale = 'log';
% xlim([1e8 1e16])
ylim([0 yl])
% ax(np).XTick = 10.^([2:2:31]);
% ax(np).XTick = 10.^([8:2:16]);
% ax(np).XTickLabel = {'10^{8}','','10^{12}','','10^{16}'}; %10.^([-2:16]);
% ax(np).TickLabelInterpreter = 'tex';
xlabel({'Supersaturation'; '(H$_2$O $wt.\%$)'},'Interpreter','Latex')
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;
hold on
xl = xlim;
plot(xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  

grid on 
set(ax(np),'XMinorTick','off')
set(ax(np),'YMinorTick','off')
ax(np).YMinorGrid = 'off';
ax(np).XMinorGrid = 'off';

% text(ax(np).XLim(1),.3,...
%     '\bf{(f)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 




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

for np = 2:(nf*nc)
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

function p2 = Conduit_plot_2(cI,cO,ax,color)

    % ylim
%     yl = 6;
%     co = get(gca,'ColorOrder');
%     color = co(2,:);

    lw = 1.25; % Line Width  
    lw2 = 1;
%     fs = 10; % Font size


    % figure('position',[50 50 1200 1200])
%     nf = 6; % Number of subplots
%     nc = 2;
%     np = 6;

    %===================================================
    % Lithostatic pressure
%     Z = linspace(0,10000,length(cO.pm))';
    plith = (cO.Z*2400*9.81+cI.pf);
  
    % Straight up pressures
    p2 = plot(ax(1),cO.pm/1e6,...
         cO.Z/1e3,...
        'LineWidth',lw,'Color',color(1,:));        
    plot(ax(1),plith/1e6,cO.Z/1e3,'--','Color',color(1,:));
    xlab = {'Pressure (MPa)'};  
    hold on
    xl = xlim(ax(1));
    plot(ax(1),xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  

    % Overpressure ratio
%     plot(ax(1),cO.pm./plith,...
%          cO.Z/1e3,...
%         'LineWidth',lw,'Color',color);        
    %===================================================
    % U

%     pp1 = plot(ax(2),cO.U,...
%                cO.Z/1e3,...
%                'LineWidth',lw,'Color',color(1,:));   

    pp2 = plot(ax(2),cO.M,...
               cO.Z/1e3,...
               'LineWidth',lw,'Color',color(1,:));    
    hold on
    xl = xlim(ax(2));
    plot(ax(2),xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  

    %===================================================
    % Decompression rate

    plot(ax(3),cO.dPdt/1e6,...
         cO.Z/1e3,...
         'LineWidth',lw,'Color',color(1,:));    
    hold on
    xl = xlim(ax(3));
    plot(ax(3),xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  
     
    %===================================================
    % Gas volume fraction

    plot(ax(4),cO.porosity * 100,...
         cO.Z/1e3,...
         'LineWidth',lw,'Color',color(1,:));  
    hold on
    xl = xlim(ax(4));
    plot(ax(4),xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  
     
    %===================================================
    % Nucleation rate
    cO.J(cO.J<1e-10) = 1e-10;
    plot(ax(5),cO.J,...
         cO.Z/1e3,...
         'LineWidth',lw,'Color',color(1,:));     
    hold on
    xl = xlim(ax(5));
    plot(ax(5),xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  

    %===================================================
    % Supersaturation
%     plot(ax(6),(cO.Cm-cO.Csol)*100,cO.Z/1e3,'LineWidth',lw,'Color',color(1,:)); 
    plot(ax(6),(cO.Cm)*100,cO.Z/1e3,'LineWidth',lw,'Color',color(1,:)); 
    plot(ax(6),(cO.Csol)*100,cO.Z/1e3,'--','LineWidth',lw,'Color',color(1,:)); 
    
%     plot(ax(6),cO.M0,...
%          cO.Z/1e3,...
%          'LineWidth',lw,'Color',color(1,:));        
    hold on
    xl = xlim(ax(6));
    plot(ax(6),xl,cO.Par.Zf*[1 1]/1e3,':','LineWidth',lw2,'Color',color(1,:));  

end