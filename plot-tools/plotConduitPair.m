function plotRunPair(d1,d2)
% plotCouplePlume(D)
% plots the detailed output of a pair of conduit runs
% d1 = struct: d1.conduitI = output of getConduitSource
%              d1.conduitO = output of Conduit_flow_with_nucleation
% Same for d2
%
% CRowell Mar 2021

    
    nrows = 1;
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
%     [~,pi] = max([d1.plumeO.hm d2.plumeO.hm]);
%     if pi==2
%         d = d1; d1 = d2; d2 = d;
%     end
    
%     if all(isfield(d1,{'conduitI','conduitO'}))
        axC = Conduit_flow_plot(d1.cI,d1.cO);
        Conduit_plot_2(d2.cI,d2.cO,axC)
        
    
end

function ax = Conduit_flow_plot(varargin)
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
%     assert(and(istable(Input),isstruct(Output)),'Input must be table, Output must be struct')
end



frag_index = find(Output.porosity>=Input.phi_frag,1);


% ylim
yl = 6;
color = [0 0 .8];
    
lw = 1.25; % Line Width  
fs = 11; % Font size
    

% figure('position',[50 50 1200 1200])
nf = 6; % Number of subplots
nc = 1;
np = 0;

%===================================================
% Pressure
np = np + 1;
ax(1) = subplot(nc,nf,np);
hold on

   
plot(Output.pm/1e6,...
     Output.Z/1e3,...
    'LineWidth',lw,'Color',color(1,:));        


% Lithostatic pressure
Z = linspace(0,10000,1000);
plot((Z*2400*9.81+Input.pf)/1e6,Z/1e3,'--','Color',color(1,:));
  


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

% text(ax(np).XLim(2)-10,.3,...
%     '\bf{(a)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 

tt = text(500*2400*9.81/1e6+10,.5,'Lithostatic','FontSize',fs);
tt.Rotation = -70; %-58;


%===================================================
% U
np = np + 1;
ax(2) = subplot(nc,nf,np);
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


% text(ax(np).XLim(1)+10,.3,...
%     '\bf{(b)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 


legend([pp1,pp2],'magma','sound','Location','southeast')

%===================================================
% Decompression rate
np = np + 1;
ax(3) = subplot(nc,nf,np);
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



% text(ax(np).XLim(2)-40,...
%     .3,...
%     '\bf{(c)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 

% grid on



%===================================================
% Gas volume fraction
np = np + 1;
ax(4) = subplot(nc,nf,np);
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
% ylabel('Depth, Z (km)')

grid on

% text(ax(np).XLim(1),.3,...
%     '\bf{(d)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 



%===================================================
% Nucleation rate
np = np + 1;
ax(5) = subplot(nc,nf,np);
hold on



Output.J(Output.J<1e-10) = 1e-10;
plot(Output.J,...
     Output.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
ax(np).XScale = 'log';
xlim([1e6 1e16])
ylim([0 yl])
ax(np).XTick = 10.^([6:4:20]);
ax(np).TickLabelInterpreter = 'tex';
xlabel({'Nucleation rate (m^{-3} s^{-1})'})
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;


% grid on

plot([1e-3 1e-3],[2.65,.8],'Color','k');



%  text(ax(np).XLim(1),.3,...
%     '\bf{(e)}','FontSize',fs,'FontWeight','bold',...
%     'VerticalAlignment', 'bottom') 


%===================================================
% BND
np = np + 1;
ax(6) = subplot(nc,nf,np);
hold on


plot(Output.M0,...
     Output.Z/1e3,...
     'LineWidth',lw,'Color',color(1,:));        


ax(np) = gca;
ax(np).YDir = 'reverse';
ax(np).XScale = 'log';
xlim([1e8 1e16])
ylim([0 yl])
ax(np).XTick = 10.^([2:2:31]);
ax(np).TickLabelInterpreter = 'tex';
xlabel({'Bubble number density (m^{-3})'})
ax(np).FontName = 'Arial';
ax(np).FontSize = fs;


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

function Conduit_plot_2(cI,cO,ax)

    % ylim
%     yl = 6;
    co = get(gca,'ColorOrder');
    color = co(2,:);

    lw = 1.25; % Line Width  
%     fs = 10; % Font size


    % figure('position',[50 50 1200 1200])
%     nf = 6; % Number of subplots
%     nc = 2;
%     np = 6;

    %===================================================
    % Pressure
%     set(ax(1),'ColorOrderIndex',2)
    plot(ax(1),cO.pm/1e6,...
     cO.Z/1e3,...
    'LineWidth',lw,'Color',color(1,:));        

    % Lithostatic pressure
    Z = linspace(0,10000,1000);
    plot(ax(1),(Z*2400*9.81+cI.pf)/1e6,Z/1e3,'--','Color',color(1,:));

    %===================================================
    % U

    pp1 = plot(ax(2),cO.U,...
               cO.Z/1e3,...
               'LineWidth',lw,'Color',color(1,:));   

    pp2 = plot(ax(2),cO.U./cO.M,...
               cO.Z/1e3,':',...
               'LineWidth',lw,'Color',color(1,:));    

    %===================================================
    % Decompression rate

    plot(ax(3),cO.dPdt/1e6,...
         cO.Z/1e3,...
         'LineWidth',lw,'Color',color(1,:));    
     
    %===================================================
    % Gas volume fraction

    plot(ax(4),cO.porosity * 100,...
         cO.Z/1e3,...
         'LineWidth',lw,'Color',color(1,:));  
     
    %===================================================
    % Nucleation rate
    cO.J(cO.J<1e-10) = 1e-10;
    plot(ax(5),cO.J,...
         cO.Z/1e3,...
         'LineWidth',lw,'Color',color(1,:));     

    %===================================================
    % BND
    plot(ax(6),cO.M0,...
         cO.Z/1e3,...
         'LineWidth',lw,'Color',color(1,:));        

end