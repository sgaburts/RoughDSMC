% Anglular scattering plots for random surface runs
close all
clear
clc


name = 'scattering.txt';

% path = '/Users/antonburtsev_1/DocumentsLocal/DARPA/new_rand_specular/rand_1/'; % texture 2
% rms = [0.1 0.5 1 2 3 5 7 9 10];
% % peak = [0.730000000000000,4.05000000000000,6.50000000000000,13.3500000000000,23.7500000000000,40.7500000000000,54.2500000000000,62.2500000000000,64.7500000000000];
% peak = [1.65000000000000,8.20000000000000,16.2000000000000,32,51,95,122,138,142]; %fwhmx

% texture 3
% path = '/Users/antonburtsev_1/DocumentsLocal/DARPA/new_rand_specular/rand_2/';
% rms = [0.1 0.5 1 2 3 5 7 9 10]; % rand_2
% % peak = [0.485000000000000,1.87500000000000,3.65000000000000,6.70000000000000,10.6500000000000,20.7500000000000,30.7500000000000,36.7500000000000,40.7500000000000]; 
% peak = [0.860000000000000,4.20000000000000,8.20000000000000,16.8000000000000,27,46,65,85,96]; % fwhmx

% texture 1
path = '/Users/antonburtsev_1/DocumentsLocal/DARPA/new_rand_specular/rand_3/';
rms = [0.1 0.5 1 2 3 5 7 9 10]; % rand_3
% peak = [1.32500000000000,6.50000000000000,12.1500000000000,29.7500000000000,43.7500000000000,63.7500000000000,69.7500000000000];
peak = [2.85000000000000,14.1000000000000,29,64,102,140,150, 158, 158]; % fwhmx

meanslope = ones(1,length(rms));
modeslope = ones(1,length(rms));



plate_x = [1];
plate_z = [1];

lines = {'-';'--';'-.';':'};
marker = ['o','d','s','^'];

% suffix = 'x magnification';
suffix = ' rms';

colormap_flag = 0;

n = 1;

colors = colororder;
newcolors = colororder; % to overwrite default color palette
% newcolors = turbo(10); % to overwrite default color palette
% newcolors = jet(10); % to overwrite default color palette

for i=1:n
%     path = ['/Users/antonburtsev/DocumentsLocal/DARPA/new_rand_specular/',num2str(cases(i)),'x_scattering.txt'];
    path = [path,name];
    
    
    T = readtable(path);
    
    theta_beam = table2array(T(2:end,1));
    Data(i).theta_beam = theta_beam;
    Fx = table2array(T(2:end,2));
    Fy = table2array(T(2:end,3));
    Fz = table2array(T(2:end,4));
    
    nn = length(theta_beam);
    
    Data(i).Drag = Fx.*sind(theta_beam) - Fy.*cosd(theta_beam);
    Data(i).A = plate_x(i)*plate_z(i) * ones(nn,1);
    Data(i).An = plate_x(i)*sind(90-theta_beam) * plate_z(i);
    
    
    Data(i).theta = table2array(T(1,5:end));
    Data(i).profile = table2array(T(2:end,5:end)); 
    

    % Normalise by maximum for each roughness case
    Data(i).profile = Data(i).profile ./ max(max(Data(i).profile));
    
%     nn = length(Data(i).profile(:,1));

    for j = 1:nn
        avg(j) = mean(Data(i).profile(j,:));
        [~,jmax] = max(Data(i).profile(j,:));
        % Correct for spurious peak at 
%         if jmax == max(theta_beam)
%             Data(i).profile(j,jmax) = 0;
%         end
        [~,jmax] = max(Data(i).profile(j,:));
        Data(i).Peak(j) = Data(i).theta(jmax);
        Data(i).diff(j) = Data(i).theta_beam(j) - Data(i).Peak(j);
        % Normalise by maximum for each roughness case and each angle
%         Data(i).profile(j,:) = Data(i).profile(j,:) ./ max(Data(i).profile(j,1:end-5));
    end
    

    %% Cartesian plot     
    
    figure(1)
    colororder(newcolors);
%     set(gca,'ColorOrderIndex',1)
    wd = 1.5;
    if i==n
        wd = 2;
    end
    for j=1:nn
        set(gca,'ColorOrderIndex',j)
        plot(Data(i).theta,Data(i).profile(j,:),cell2mat(lines(i)),'LineWidth',wd); hold on
        

%         data = Data(i).profile(j,:)./(2*pi*sin(deg2rad(Data(i).theta)) * deg2rad(1));
%         data = data./max(data);
%         plot(Data(i).theta,data,cell2mat(lines(i)),'LineWidth',wd); hold on
    end
    
    plot(Data(i).theta,cosd(Data(i).theta),'k-.','LineWidth',3); hold on
    
    xlabel({'$\Theta ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    ylabel({'Normalized Flux'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
%     legend show

    if colormap_flag==1
    colormap(turbo);
    c = colorbar('FontSize',18);
    c.Location = 'eastoutside';
    c.Label.String = '$\Theta_b$';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    c.LineWidth = 1.5;
    c.Label.Rotation = 0;
    c.Label.FontSize = 20;
    caxis([min(theta_beam), max(theta_beam)])
    set(gca,'xTick',[0:10:90])
    end
    
    %% Polar plot
    for j = 1:nn
        % Normalise by maximum for each roughness case and each angle
        Data(i).profile(j,:) = Data(i).profile(j,:) ./ max(Data(i).profile(j,1:end-5));
    end
    
%     newcolors = turbo(10);
    
    figure(2)
    colororder(newcolors);
%     set(gca,'ColorOrderIndex',1)
    
    for j=1:nn
        
        % Smooth and renormalise
        Data(i).profile(j,:) = smoothdata(Data(i).profile(j,:),"gaussian",6);
        Data(i).profile(j,:) = Data(i).profile(j,:) ./ max(max(Data(i).profile(j,:)));
        
        set(gca,'ColorOrderIndex',j)
        polarplot(deg2rad(Data(i).theta),Data(i).profile(j,:),cell2mat(lines(i)),'LineWidth',2); hold on

%         data = Data(i).profile(j,:)./(2*pi*sin(deg2rad(Data(i).theta)) * deg2rad(1));
%         data = data./max(data);
%         polarplot(deg2rad(Data(i).theta),data,cell2mat(lines(i)),'LineWidth',2); hold on
    end
    
    polarplot(deg2rad(Data(i).theta),cosd(Data(i).theta),'k-.','LineWidth',3); hold on
    
%     ax = gca;
%     box on
%     ax.ThetaLim = [0 90];
%     ax.ThetaDir = 'clockwise';
%     ax.ThetaZeroLocation = "top";
%     ax.ThetaTick = [0:10:90];
%     ax.ThetaMinorTick = 'on';
%     box on;
%     set(gca,'FontSize',18)      % Axis and labels fontsise
%     set(gca,'linewidth',1.5)    % Axis line width4
%     set(gca,'TickLabelInterpreter','latex')
%     ax.RLim = [0 1];

    ax = gca;
    box on
    ax.ThetaLim = [-90 90];
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = "top";
    ax.ThetaTick = [-90:10:90];
    ax.ThetaMinorTick = 'on';
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    
    if colormap_flag==1
    colormap(turbo);
    c = colorbar('FontSize',18);
    c.Location = 'eastoutside';
    c.Label.String = '$\Theta_b$';
    c.Label.Interpreter = 'latex';
    c.TickLabelInterpreter = 'latex';
    c.LineWidth = 1.5;
    c.Label.Rotation = 0;
    c.Label.FontSize = 20;
    caxis([min(theta_beam), max(theta_beam)])
    end
    
    %%
    
    j = i;
%     j = 3;
    
%     figure(3)
%     colororder(colors);
%     set(gca,'ColorOrderIndex',j)
%     plot(rms,Data(i).Drag ./ Data(i).An,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8); hold on
%     
% %     set(gca,'ColorOrderIndex',j)
% %     plot(rms,Data(i).Drag(2:end-1) ./ Data(i).An(2:end-1),[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8); hold on
% %     set(gca,'ColorOrderIndex',j)
% %     plot([rms(1) rms(end)],[Data(i).Drag(1) ./ Data(i).An(1) Data(i).Drag(1) ./ Data(i).An(1)],'-','LineWidth',2,'MarkerSize',8);
% %     set(gca,'ColorOrderIndex',j)
% %     plot([rms(1) rms(end)],[Data(i).Drag(end) ./ Data(i).An(end) Data(i).Drag(end) ./ Data(i).An(end)],'--','LineWidth',2,'MarkerSize',8);
%     
%     xlabel({'rms'},'Interpreter','latex','FontSize',20)
%     ylabel({'$D/A_n ~(N/m^2)$'},'Interpreter','latex','FontSize',20)
%     box on;
%     set(gca,'FontSize',18)      % Axis and labels fontsise
%     set(gca,'linewidth',1.5)    % Axis line width4
%     set(gca,'TickLabelInterpreter','latex')
%     
    figure(4)
    colororder(colors);
%     set(gca,'ColorOrderIndex',j)
    plot(rms,Data(i).Drag ./ Data(i).A,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8); hold on
%     set(gca,'ColorOrderIndex',j)
%     plot([rms(1) rms(end)]*1000,[Data(i).Drag(1) ./ Data(i).A(1) Data(i).Drag(1) ./ Data(i).A(1)],'-','LineWidth',2,'MarkerSize',8);
%     set(gca,'ColorOrderIndex',j)
%     plot([rms(1) rms(end)]*1000,[Data(i).Drag(end) ./ Data(i).A(end) Data(i).Drag(end) ./ Data(i).A(end)],'--','LineWidth',2,'MarkerSize',8);
    
    % xlabel({'RMS ($n m$)'},'Interpreter','latex','FontSize',20)
    xlabel({'$S_q$ ($n m$)'},'Interpreter','latex','FontSize',20)
    ylabel({'$D/A ~ (N/m^2)$'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    


%     
    
    %
    figure(7)
    colororder(colors);
    plot(peak,Data(i).Drag ./ Data(i).A,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8); hold on
    
    xlabel({'peak $\theta_y$'},'Interpreter','latex','FontSize',20)
    ylabel({'$D/A ~ (N/m^2)$'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')


end