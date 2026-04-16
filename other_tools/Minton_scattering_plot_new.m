close all
clear
clc

% Plot Minton scattering data
path = '/Users/antonburtsev_1/DocumentsLocal/DARPA/UCB_Scattering_1/';

% cases = {'UTAustin','UTDallas','GTRI'};
cases = {'UTDallas','UTAustin'};

theta_beam = [70 50 30]; % all angles
angles = [1 2 3];   % angles to plot

markers = ['o','s','d'];

Tw = 300;
kb = 1.380649*10^(-23); %J/K
Na = 6.02214076*10^23;

E_in = 495.84234*1000; % for O
% E_in = 996*1000; % for O2


style = 2;   % 1 (Color-angle, marker-case); 2 (Color-case, marker-angle); 3 change both;
normalise = 1;

leg_font = 20; %18


colors = colororder;
n1 = length(cases);
n2 = length(theta_beam);
% colors(n1+1:n1+n2,:) = 0;

for caseind = 1:n1
%     for angleind = 1:n2
        for angleind = angles(1):angles(end)
        
        if style == 1
            % Color - angle, marker - case
            icol  = angleind;
            imark = caseind;
        elseif style == 2
            % Color - case, marker - angle
            icol  = caseind;
            imark = angleind;
        elseif style == 3
            % Color - angle, marker - angle
            if length(angles) == 1
            icol  = caseind;
            imark = caseind;
            else
            icol  = angleind;
            imark = angleind;
            end
        end
            

    T = readtable([path,char(cases(caseind)),'.xlsx'],'Sheet',angleind);
    A = table2array(T);
    theta = A(:,1); 
    E_out = A(:,2); % 2-O, 3-O2
    count = A(:,4); % 4-O, 5-O2
    
    if normalise == 1
        count = count./max(count);
    end
        
    E_w = 3/2 * Na * kb * Tw;
    E_out = E_out*1000;
    alpha = (E_in - E_out)./(E_in - E_w);
   
    xi = 180 - (theta_beam(angleind) + theta);
    label = [char(cases(caseind)),' $\Theta_b = $',num2str(theta_beam(angleind))];
    
%     figure(1)
%     plot(theta,count); hold on
    
    figure(2)
    scatter(theta,E_out./E_in,100,markers(imark),'LineWidth',2,'MarkerEdgeColor',colors(icol,:),'MarkerFaceColor',colors(icol,:),'MarkerFaceAlpha',0.2,'MarkerFaceAlpha',0.2,'DisplayName',label); hold on
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$\Theta$ (deg)','Interpreter','latex')
    ylabel('$\langle E_f \rangle/\langle E_i \rangle$','Interpreter','latex')
    legend('show','Interpreter','latex','box','off','FontSize',leg_font)
    xlim([-20 90])
    ylim([0 1])
    marker_size_fix(10);
    
    figure(3)
    box on
    polarscatter(deg2rad(theta),count,100,markers(imark),'LineWidth',2,'MarkerEdgeColor',colors(icol,:),'MarkerFaceColor',colors(icol,:),'MarkerFaceAlpha',0.2,'MarkerFaceAlpha',0.2,'DisplayName',label); hold on
    ax = gca;
    ax.ThetaLim = [0 90];
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = "top";
    ax.ThetaTick = [0:10:90];
    ax.ThetaMinorTick = 'on';
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    if normalise == 1
        ax.RLim = [0 1.05];
    end
    legend('show','Interpreter','latex','box','off','FontSize',leg_font)
    marker_size_fix(10);
    
    figure(4)
    scatter(theta,alpha,100,markers(imark),'LineWidth',2,'MarkerEdgeColor',colors(icol,:),'MarkerFaceColor',colors(icol,:),'MarkerFaceAlpha',0.2,'MarkerFaceAlpha',0.2,'DisplayName',label); hold on
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$\Theta$ (deg)','Interpreter','latex')
    ylabel('$\alpha$','Interpreter','latex')
    legend('show','Interpreter','latex','box','off','FontSize',leg_font)
    xlim([-20 90])
    ylim([0 1])
    marker_size_fix(10);
    
    figure(5)
    scatter(xi,(E_in-E_out)./E_in,100,markers(imark),'LineWidth',2,'MarkerEdgeColor',colors(icol,:),'MarkerFaceColor',colors(icol,:),'MarkerFaceAlpha',0.2,'DisplayName',label); hold on
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$\chi$ (deg)','Interpreter','latex')
    ylabel('$(\langle E_i \rangle - \langle E_f \rangle)/\langle E_i \rangle$','Interpreter','latex')
    legend('show','Interpreter','latex','box','off','FontSize',leg_font)
    xlim([20 140])
    ylim([0 1])
    marker_size_fix(10);
    
    end
    
end


% figure(4)
% hLegend = findobj(gcf, 'Type', 'Legend'); %find existing legend
% for i=1:5
% hLegendEntry = hLegend.EntryContainer.NodeChildren(i);
% hLegendIconLine = hLegendEntry.Icon.Transform.Children.Children;
% hLegendIconLine(end).Size=10;
% end

%%
function [] = marker_size_fix(sz)
%marker_size_fix Change size of scatter markers in the legend
%   Detailed explanation goes here
hLegend = findobj(gcf, 'Type', 'Legend'); %find existing legend
n = length(hLegend.String); %number of legend entries
for i=1:n
    hLegendEntry = hLegend.EntryContainer.NodeChildren(i);
    hLegendIconLine = hLegendEntry.Icon.Transform.Children.Children;
    hLegendIconLine(end).Size = sz;
end
end

