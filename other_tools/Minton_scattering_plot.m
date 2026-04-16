close all
clear
clc

% Plot Minton scattering data
path = '/Users/antonburtsev_1/Documents/University/AUSTIN/DARPA/Minton_paper_data/';
% cases = {'Austin_Scattering_1_70_incedent_angle.csv','Austin_Scattering_1_50_incedent_angle.csv','Austin_Scattering_1_30_incedent_angle.csv'};
cases = {'Saphire-70deg.csv','BFS30-70deg.csv','Coated_unexposed_POSS.csv','Uncoated_exposed_POSS.csv','Austin_Scattering_1_70_incedent_angle.csv'};
% cases = {'Austin_Scattering_1_70_incedent_angle.csv'};


Tw = 300;
kb = 1.380649*10^(-23); %J/K
Na = 6.02214076*10^23;
% kb = 8.617333262*10^(-5); %eV/K

% E_in = 495.84234*1000; % for O
% E_in = 996*1000; % for O2

theta_beam = [70 70 70 70 70];
% markers = ['o','s','d','^','d'];
markers = ['v','^','d','s','o'];
labels = {'Sapphire','BFS-30 (AO exposed)','$\mathrm{Al_2 O_3}$ coated POSS','POSS (AO exposed)','Ge coated Si'};

n1 = length(cases);
colors = colororder;
colors(3,:) = [];
colors = colors(1:n1,:);
colors = flipud(colors);
colors(end-1:end,:) = flipud(colors(end-1:end,:));

for i=1:n1

    C = readcell([path,char(cases(i))]);
    E_in = C{1,1}; % <- read E_i from header
    A = importdata([path,char(cases(i))],',',2);
    theta = A.data(:,1);
    E_out = A.data(:,2);
    if max(max(E_out)<10) % <- ratio of E_f/E_i
        E_out = E_out*E_in;
    end
    E_out = E_out/1000;
    count = A.data(:,4); % 4-O, 5-O2
    
    E_w = 3/2 * Na * kb * Tw;
    E_out = E_out*1000;
    alpha = (E_in - E_out)./(E_in - E_w);
   
    xi = 180 - (theta_beam(i) + theta);

    
    figure(1)
    scatter(theta,count./max(count),100,markers(i),'LineWidth',2,'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerFaceAlpha',0.2,'MarkerFaceAlpha',0.2,'DisplayName',char(labels(i))); hold on
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$\Theta$ (deg)','Interpreter','latex')
    ylabel('Normalised Flux','Interpreter','latex')
    lgd = legend('show','Interpreter','latex','box','off','FontSize',18);
    ylim([0 1.05])
    marker_size_fix(lgd,10);
    
    

    figure(2)
%     plot(theta,E_out./E_in,markers(i),'LineWidth',2,'MarkerSize',10); hold on
    scatter(theta,E_out./E_in,100,markers(i),'LineWidth',2,'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerFaceAlpha',0.2,'MarkerFaceAlpha',0.2,'DisplayName',char(labels(i))); hold on
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$\Theta$ (deg)','Interpreter','latex')
    ylabel('$\langle E_f \rangle/\langle E_i \rangle$','Interpreter','latex')
    lgd = legend('show','Interpreter','latex','box','off','FontSize',18);
    marker_size_fix(lgd,10);
    
    figure(3)
    box on
%     polarplot(deg2rad(theta),count./max(count),markers(i),'LineWidth',2,'MarkerSize',10); hold on
    polarscatter(deg2rad(theta),count./max(count),150,markers(i),'LineWidth',2,'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerFaceAlpha',0.2,'MarkerFaceAlpha',0.2,'DisplayName',char(labels(i))); hold on
    ax = gca;
    ax.ThetaLim = [0 90];
    % ax.ThetaLim = [-90 90];
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = "top";
    ax.ThetaTick = [-90:10:90];
    ax.ThetaMinorTick = 'on';
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    ax.RLim = [0 1.05];
    lgd = legend('show','Interpreter','latex','box','off','FontSize',18);
    marker_size_fix(lgd,10);
    
    figure(4)
%     plot(theta,alpha,markers(i),'LineWidth',2,'MarkerSize',10); hold on
    scatter(theta,alpha,100,markers(i),'LineWidth',2,'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerFaceAlpha',0.2,'MarkerFaceAlpha',0.2,'DisplayName',char(labels(i))); hold on
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$\Theta$ (deg)','Interpreter','latex')
    ylabel('$\alpha$','Interpreter','latex')
    lgd = legend('show','Interpreter','latex','box','off','FontSize',18);
    marker_size_fix(lgd,10);
    
    figure(5)
    scatter(xi,(E_in-E_out)./E_in,100,markers(i),'LineWidth',2,'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',colors(i,:),'MarkerFaceAlpha',0.2,'DisplayName',char(labels(i))); hold on
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$\chi$ (deg)','Interpreter','latex')
    ylabel('$(\langle E_i \rangle - \langle E_f \rangle)/\langle E_i \rangle$','Interpreter','latex')
    lgd = legend('show','Interpreter','latex','box','off','FontSize',18);
    marker_size_fix(lgd,10);
    
end





%%
function[] =  marker_size_fix(hLegend, sz)
drawnow;
n = numel(hLegend.String);
for i = 1:n
    hLegendEntry = hLegend.EntryContainer.NodeChildren(i);
    hLegendIconLine = hLegendEntry.Icon.Transform.Children.Children;
    hLegendIconLine(end).Size = sz;
end

end


