close all
clear
clc

resultspath = '/Users/antonburtsev_1/DocumentsLocal/DARPA/AFM_drag/';


% resultspath = '/Users/antonburtsev_1/DocumentsLocal/DARPA/';
% 
% resultspath = '/Users/antonburtsev_1/DocumentsLocal/DARPA/new_rand_specular/rand_2/';

% cases = {'Kapton','Kapton_1h_Exp','Ge_Si','Ge_Si_LSCM','Ge_Si_Exp_LSCM'};
% cases = {'Kapton','Kapton_1h_Exp','Ge_Si','Ge_Si_1h_Exp','Kapton_1h_Exp_2'};
% cases = {'Ge_Si_LSCM','Ge_Si_LSCM_2','Ge_Si_Exp_LSCM','Ge_Si_Exp_LSCM_2'};
% cases = {'new/Ge_Si_AFM','Ge_Si_LSCM','new/Ge_Si_AFM_EXP_UCB','Ge_Si_Exp_LSCM'};
% cases = {'new/Ge_Si_AFM','Ge_Si_LSCM'};
cases = {'new/Ge_Si_AFM_newcll70','new/Ge_Si_AFM'};


% cases = {'2.5x_scattering','5x_scattering','10x_scattering','20x_scattering'};
% cases = {'thin_specular_2','thin_diffuse_2','thin_specular_3','thin_diffuse_3'};


% plate_x = 2*[0.0001762 0.0001762 0.0001762 0.0001762 0.0001762 0.0001762];
% plate_z = 2*[0.0001320 0.0001320 0.0001320 0.0001320 0.0001320 0.0001320];

% For AFM
plate_x = [4.9805 4.9805 4.9805 4.9805 4.9805 4.9805 4.9805 4.9805 4.9805 4.9805]*1e-6;
plate_z = [4.9805 4.9805 4.9805 4.9805 4.9805 4.9805 4.9805 4.9805 4.9805 4.9805]*1e-6;

% plate_x = [1, 1, 1, 1];
% plate_z = [1, 1, 1, 1];

% plate_x = [4.9805*1e-6 4.9805*1e-6 4.9805*1e-6 2*0.0001762 2*0.0001762];
% plate_z = [4.9805*1e-6 4.9805*1e-6 4.9805*1e-6 2*0.0001320 2*0.0001320];


exp_path = '/Users/antonburtsev_1/Documents/University/AUSTIN/DARPA/Minton_paper_data/';
exp_cases = {'Austin_Scattering_1_70_incedent_angle.csv','Austin_Scattering_1_50_incedent_angle.csv'};

lines = {'-';'-.';'--';'-.';'-.';'-.'};
marker = ['o','d','s','^','v','>'];

suffix = 'x magnification';

colormap_flag = 0;

n = length(cases);

colors = colororder;
% newcolors = turbo(10); % to overwrite default color palette
newcolors = colors;

m = 2.65e-26;

m = 1/(0.74454/(2.65e-26) + 0.2423/(4.65e-26) + 0.01315/(5.31e-26));

nrho = 1.45e15;
U = 7896.14402;

plot_angle = 2; % 1 beam angle, 2 AoA

%% Load smooth cases

smooth_plate_x = 2*1.37882;
smooth_plate_z = 2*1.03411;
% smooth_plate_x = 1;
% smooth_plate_z = 1;
A = smooth_plate_x*smooth_plate_z;

% path = ['/Users/antonburtsev_1/DocumentsLocal/DARPA/thin_specular.txt'];
path = ['/Users/antonburtsev_1/DocumentsLocal/DARPA/results/smooth_specular_scattering.txt'];
T = readtable(path);
theta_beam_smooth = table2array(T(2:end,1));
Fx = table2array(T(2:end,2));
Fy = table2array(T(2:end,3));
Fz = table2array(T(2:end,4));
An = smooth_plate_x*sind(90-theta_beam_smooth) * smooth_plate_z;
Drag_spec = Fx.*sind(theta_beam_smooth) - Fy.*cosd(theta_beam_smooth);
Lift_spec = -Fx.*cosd(theta_beam_smooth) - Fy.*sind(theta_beam_smooth);


% path = ['/Users/antonburtsev_1/DocumentsLocal/DARPA/thin_diffuse.txt'];
path = ['/Users/antonburtsev_1/DocumentsLocal/DARPA/results/smooth_diffuse_scattering.txt'];
T = readtable(path);
theta_beam_smooth = table2array(T(2:end,1));
Fx = table2array(T(2:end,2));
Fy = table2array(T(2:end,3));
Fz = table2array(T(2:end,4));
An = smooth_plate_x*sind(90-theta_beam_smooth) * smooth_plate_z;
Drag_diff = Fx.*sind(theta_beam_smooth) - Fy.*cosd(theta_beam_smooth);
Lift_diff = -Fx.*cosd(theta_beam_smooth) - Fy.*sind(theta_beam_smooth);

% path = ['/Users/antonburtsev/DocumentsLocal/DARPA/results/smooth_cll_scattering.txt'];
% smooth_plate_x = 0.0970899;
% smooth_plate_z = 0.0725795;
% A2 = smooth_plate_x*smooth_plate_z;
% 
% T = readtable(path);
% theta_beam_smooth_2 = table2array(T(2:end,1));
% Fx = table2array(T(2:end,2));
% Fy = table2array(T(2:end,3));
% Fz = table2array(T(2:end,4));
% An = smooth_plate_x*sind(90-theta_beam_smooth_2) * smooth_plate_z;
% Drag_cll = Fx.*sind(theta_beam_smooth_2) - Fy.*cosd(theta_beam_smooth_2);

%% Load DM data

MD = readmatrix('/Users/antonburtsev_1/Documents/University/AUSTIN/DARPA/md_lit_data.xlsx');

%% Load EXP data

n_exp = length(exp_cases);

for i=1:n_exp
    C = readcell([exp_path,char(exp_cases(i))]);
    exp_E_in = C{1,1}; % <- read E_i from header
    B = importdata([exp_path,char(exp_cases(i))],',',2);
    exp(i).theta = B.data(:,1);
    exp(i).E_out = B.data(:,2);
    if max(max(exp(i).E_out)<10) % <- ratio of E_f/E_i
        exp(i).E_out = exp(i).E_out*E_in;
    end
    exp(i).E_out = exp(i).E_out/1000;
    exp(i).count = B.data(:,4); % 4-O, 5-O2
    exp(i).count = exp(i).count./max(exp(i).count)
end


%% Plotting
for i=1:n
    % path = [resultspath,char(cases(i)),'_scattering.txt'];
    path = [resultspath,char(cases(i)),'.txt'];
    
    T = readtable(path);
    
    theta_beam = table2array(T(2:end,1));
    Data(i).theta_beam = theta_beam;
    Fx = table2array(T(2:end,2));
    Fy = table2array(T(2:end,3));
    Fz = table2array(T(2:end,4));
    
    Data(i).Drag = Fx.*sind(theta_beam) - Fy.*cosd(theta_beam);
    Data(i).Lift = -Fx.*cosd(theta_beam) - Fy.*sind(theta_beam);
%     Data(i).Drag = Fy;

    Data(i).A = plate_x(i)*plate_z(i);
    Data(i).An = plate_x(i)*sind(90-theta_beam) * plate_z(i);
    
    
    Data(i).theta = table2array(T(1,5:end));
    Data(i).profile = table2array(T(2:end,5:end)); 
    

    % Normalise by maximum for each roughness case
    Data(i).profile = Data(i).profile ./ max(max(Data(i).profile));

    for j = 1:length(Data(i).profile(:,1))
        [~,jmax] = max(Data(i).profile(j,:));
        % Correct for spurious peak at 
        if jmax == max(theta_beam)
            Data(i).profile(j,jmax) = 0;
        end
        [~,jmax] = max(Data(i).profile(j,:));
        Data(i).Peak(j) = Data(i).theta(jmax);
        Data(i).diff(j) = Data(i).theta_beam(j) - Data(i).Peak(j);
        % Normalise by maximum for each roughness case and each angle
%         Data(i).profile(j,:) = Data(i).profile(j,:) ./ max(Data(i).profile(j,1:end-5));
    end
    
    % nn = length(Data(i).Drag);
    [tf, idx] = ismember(Data(i).theta_beam, theta_beam_smooth);

    %% Cartesian plot 

    for j = 1:length(Data(i).profile(:,1))
        % Normalise by maximum for each roughness case and each angle
        Data(i).profile(j,:) = Data(i).profile(j,:) ./ max(Data(i).profile(j,1:end-5));
    end
    
    figure(1)
    colororder(newcolors);
    set(gca,'ColorOrderIndex',1)
    wd = 1.5;
    if i==n
        wd = 2;
    end
    plot(Data(i).theta,Data(i).profile,cell2mat(lines(i)),'LineWidth',wd, 'DisplayName',char(cases(i))); hold on
    
    for j = 1:length(Data(i).profile(:,1))
        th_max(j) = Data(i).theta(find(Data(i).profile(j,:) == max(Data(i).profile(j,:))));
    end

    if i==n
        for j=1:n_exp
            plot(exp(j).theta,exp(j).count,marker(j),'LineWidth',3,'Color','k','MarkerSize',8,'DisplayName',char(exp_cases(j))); hold on
        end
    end
    
    xlabel({'$\Theta ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    ylabel({'Normalized Flux'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
%     legend show
    xlim([0 90]); set(gca,'xTick',[0:10:90])

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
    for j = 1:length(Data(i).profile(:,1))
        % Normalise by maximum for each roughness case and each angle
        Data(i).profile(j,:) = Data(i).profile(j,:) ./ max(Data(i).profile(j,1:end-5));
    end

    
%     newcolors = turbo(10);
    
    figure(2)
    colororder(newcolors);
    set(gca,'ColorOrderIndex',1)
    box on
    polarplot(deg2rad(Data(i).theta),Data(i).profile,cell2mat(lines(i)),'LineWidth',3, 'DisplayName',char(cases(i))); hold on

    if i==n
        for j=1:n_exp
            polarplot(deg2rad(exp(j).theta),exp(j).count,marker(j),'LineWidth',3,'Color','k','MarkerSize',8,'DisplayName',char(exp_cases(j))); hold on
        end
    end

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
    ax.RLim = [0 1];
    
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
    
    figure(3)
    colororder(colors);
    set(gca,'ColorOrderIndex',i)
    plot(Data(i).theta_beam,Data(i).Drag ./ Data(i).An,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8); hold on
    
    xlabel({'$\Theta_b ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    ylabel({'$D/A_n ~(N/m^2)$'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlim([0 90]); set(gca,'xTick',[0:10:90])

    
    
    figure(4)
    colororder(colors);
    if i==1
%         plot(theta_beam_smooth,Drag_spec./A,'k-o','LineWidth',3,'MarkerSize',8, 'DisplayName','Smooth Specular'); hold on
%         plot(theta_beam_smooth,Drag_diff./A,'k--o','LineWidth',3,'MarkerSize',8, 'DisplayName','Smooth Diffuse'); hold on
        plot(theta_beam_smooth,Drag_spec./A,'k-o','LineWidth',3,'MarkerSize',8,'HandleVisibility','off'); hold on
        plot(theta_beam_smooth,Drag_diff./A,'k--o','LineWidth',3,'MarkerSize',8,'HandleVisibility','off'); hold on
        
%         plot(theta_beam_smooth_2,Drag_cll./A2,'k--o','LineWidth',3,'MarkerSize',8,'HandleVisibility','off'); hold on
    end
    set(gca,'ColorOrderIndex',i)
    plot(Data(i).theta_beam,Data(i).Drag ./ Data(i).A,[marker(i),cell2mat(lines(i))],'LineWidth',3,'MarkerSize',10, 'DisplayName',char(cases(i))); hold on

    Data(i).Drag ./ Data(i).A
    
    xlabel({'$\Theta_b ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    ylabel({'$D/A ~ (N/m^2)$'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    xlim([0 90]); set(gca,'xTick',[0:10:90])
    
    leg = legend;
    leg.Interpreter = 'Latex';
    
    xlim([45 90])
    
%     spec_diff(:,i) = (Data(i).Drag ./ Data(i).A) - (Drag_spec./A);
    
%     figure(5)
%     colororder(colors);
%     set(gca,'ColorOrderIndex',i)
%     plot(Data(i).theta_beam, Data(i).diff,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8,...
%         'DisplayName',[num2str(cases(i)),suffix]); hold on
%     xlabel({'$\Theta_b ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
%     ylabel({'$\Theta_{spec} - \Theta ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
%     box on;
%     set(gca,'FontSize',18)      % Axis and labels fontsise
%     set(gca,'linewidth',1.5)    % Axis line width4
%     set(gca,'TickLabelInterpreter','latex')
%     legend('Interpreter','latex','FontSize',18,'Location','northeast','NumColumns',1,'box','off');
%     xlim([0 90]); set(gca,'xTick',[0:10:90])

%     % Difference form specular or cll
%     figure(6)
%     set(gca,'ColorOrderIndex',i)
%     plot_data = (Data(i).Drag ./ Data(i).A) - (Drag_spec(6:end-1)./A);    %from specular smooth
% %     plot_data = abs((Data(i).Drag ./ Data(i).A) - (Drag_cll./A2));                %from cll smooth
%     plot(Data(i).theta_beam,plot_data,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8, 'DisplayName',char(cases(i))); hold on
% 
%     
%     xlabel({'$\Theta_b ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
%     ylabel({'$(D - D_{\textrm{smooth}})/A ~ (N/m^2)$'},'Interpreter','latex','FontSize',20)
%     set(gca,'YScale','log')
%     box on;
%     set(gca,'FontSize',18)      % Axis and labels fontsise
%     set(gca,'linewidth',1.5)    % Axis line width4
%     set(gca,'TickLabelInterpreter','latex')
%     xlim([0 90]); set(gca,'xTick',[0:10:90])
%     
%     leg = legend;
%     leg.Interpreter = 'Latex';
%     
%     xlim([45 90])
   

    % figure(6)
    % A = smooth_plate_x * smooth_plate_z;
    % An = smooth_plate_x*sind(90-theta_beam_smooth) * smooth_plate_z;
    % % spec_diff(:,i) = (Data(i).Drag ./ Data(i).An) - (Drag_spec(end-nn:end-1)./An(end-nn:end-1));
    % spec_diff(:,i) = (Data(i).Drag ./ Data(i).An) - (Drag_spec(idx)./An(idx));
    % 
    % 
    % colororder(colors);
    % set(gca,'ColorOrderIndex',i)
    % plot(Data(i).diff,spec_diff(:,i),[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8); hold on
    % xlabel({'$\Theta_{spec} - \Theta_{max} ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    % ylabel({'$(D - D_{spec})/An $'},'Interpreter','latex','FontSize',20)
    % box on;
    % set(gca,'FontSize',18)      % Axis and labels fontsise
    % set(gca,'linewidth',1.5)    % Axis line width4
    % set(gca,'TickLabelInterpreter','latex')
    % legend('Interpreter','latex','FontSize',18,'Location','northeast','NumColumns',1,'box','off');
    % 

    % Maxwell theory
    alphas = [0:1:90];
    [Cl_md,Cd_md] = maxwell_plate(alphas,1,U,2.65e-26,300,300); % Maxwell diffuse
    [Cl_ms,Cd_ms] = maxwell_plate(alphas,0,U,2.65e-26,300,300); % Maxwell specular


    if plot_angle == 1
        angle = Data(i).theta_beam;
        angle_smooth = theta_beam_smooth; % Beam angle
        MD_angle = MD(:,1);
        alphas = 90 - alphas;
    elseif plot_angle == 2
        angle = 90-Data(i).theta_beam;
        angle_smooth = 90-theta_beam_smooth; % AoA
        MD_angle = 90-MD(:,1);
    end
    


    % Cd plot
    figure(10)
    Cd = Data(i).Drag ./(0.5*m*nrho*U^2*Data(i).A);
    Cd_spec = Drag_spec./(0.5*m*nrho*U^2*A);
    Cd_diff = Drag_diff./(0.5*m*nrho*U^2*A);
    colororder(colors);
    set(gca,'ColorOrderIndex',i)
    plot(angle,Cd,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8,'DisplayName','Ge Si Exposed DSMC'); hold on
    if i==1
        plot(angle_smooth,Cd_diff,'r--','LineWidth',2,'MarkerSize',8,'DisplayName','Smooth plate diffuse DSMC'); hold on
        plot(angle_smooth,Cd_spec,'r-','LineWidth',2,'MarkerSize',8,'DisplayName','Smooth plate specular DSMC'); hold on
        plot(MD_angle,MD(:,2),'ko','LineWidth',2,'MarkerSize',8,'DisplayName','MD amorphous $\mathrm{Al_2 O_3}$'); hold on
        plot(alphas,Cd_md,'k--','LineWidth',2,'MarkerSize',8,'DisplayName','Maxwell model diffuse'); hold on
        plot(alphas,Cd_ms,'k-','LineWidth',2,'MarkerSize',8,'DisplayName','Maxwell model specular'); hold on
    end
    if plot_angle == 1
        xlabel({'$\Theta_b ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    elseif plot_angle == 2
        xlabel({'$\alpha ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    end
    ylabel({'$C_d$'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    legend('Interpreter','latex','FontSize',18,'Location','northeast','NumColumns',1,'box','off');
    xlim([0 90])


    % Cl plot
    figure(11)
    Cl = Data(i).Lift ./(0.5*m*nrho*U^2*Data(i).A);
    Cl_spec = Lift_spec./(0.5*m*nrho*U^2*A);
    Cl_diff = Lift_diff./(0.5*m*nrho*U^2*A);
    colororder(colors);
    set(gca,'ColorOrderIndex',i)
    plot(angle,Cl,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8,'DisplayName','Ge Si Exposed DSMC'); hold on
    if i==1
        plot(angle_smooth,Cl_diff,'r--','LineWidth',2,'MarkerSize',8,'DisplayName','Smooth plate diffuse DSMC'); hold on
        plot(angle_smooth,Cl_spec,'r-','LineWidth',2,'MarkerSize',8,'DisplayName','Smooth plate specular DSMC'); hold on
        plot(MD_angle,MD(:,3),'ko','LineWidth',2,'MarkerSize',8,'DisplayName','MD amorphous $\mathrm{Al_2 O_3}$'); hold on
        plot(alphas,Cl_md,'k--','LineWidth',2,'MarkerSize',8,'DisplayName','Maxwell model diffuse'); hold on
        plot(alphas,Cl_ms,'k-','LineWidth',2,'MarkerSize',8,'DisplayName','Maxwell model specular'); hold on
    end
    if plot_angle == 1
        xlabel({'$\Theta_b ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    elseif plot_angle == 2
        xlabel({'$\alpha ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    end
    ylabel({'$C_l$'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    legend('Interpreter','latex','FontSize',18,'Location','northeast','NumColumns',1,'box','off');
    xlim([0 90])

    % Cl/Cd plot
    figure(12)
    LD = Cl./Cd;
    LD_spec = Cl_spec./Cd_spec;
    LD_diff = Cl_diff./Cd_diff;
    colororder(colors);
    set(gca,'ColorOrderIndex',i)
    plot(angle,LD,[marker(i),cell2mat(lines(i))],'LineWidth',2,'MarkerSize',8,'DisplayName','Ge Si Exposed DSMC'); hold on
    if i==1
        plot(angle_smooth,LD_diff,'r--','LineWidth',2,'MarkerSize',8,'DisplayName','Smooth plate diffuse DSMC'); hold on
        plot(angle_smooth,LD_spec,'r-','LineWidth',2,'MarkerSize',8,'DisplayName','Smooth plate specular DSMC'); hold on
        plot(MD_angle,MD(:,3)./MD(:,2),'ko','LineWidth',2,'MarkerSize',8,'DisplayName','MD amorphous $\mathrm{Al_2 O_3}$'); hold on
        plot(alphas,Cl_md./Cd_md,'k--','LineWidth',2,'MarkerSize',8,'DisplayName','Maxwell model diffuse'); hold on
        plot(alphas,Cl_ms./Cd_ms,'k-','LineWidth',2,'MarkerSize',8,'DisplayName','Maxwell model specular'); hold on
    end
    if plot_angle == 1
        xlabel({'$\Theta_b ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    elseif plot_angle == 2
        xlabel({'$\alpha ~ (^\circ)$'},'Interpreter','latex','FontSize',20)
    end
    ylabel({'$C_l/C_d$'},'Interpreter','latex','FontSize',20)
    box on;
    set(gca,'FontSize',18)      % Axis and labels fontsise
    set(gca,'linewidth',1.5)    % Axis line width4
    set(gca,'TickLabelInterpreter','latex')
    legend('Interpreter','latex','FontSize',18,'Location','northeast','NumColumns',1,'box','off');
    ylim([0 2])
    xlim([0 90])


end

%%
function [Cl,Cd] = maxwell_plate(alpha,sigma,u,m,T,Tw)

kb = 1.380649e-23;
a = deg2rad(alpha);
s = u./sqrt((2*kb*T)/m);
e = 1-sigma;

Cl = 4*e./(s.*sqrt(pi)) .* sin(a).*cos(a) ...
     .* exp(-s.^2 .* sin(a).^2) ...
   + cos(a)./s.^2 .* (1 + e.*(1 + 4*s.^2 .* sin(a).^2)) ...
     .* erf(s.*sin(a)) ...
   + (1-e)./s .* sqrt(pi) .* sin(a).*cos(a) ...
     .* sqrt(Tw/T);

Cd = 2*(1 - e.*cos(2*a))./(s.*sqrt(pi)) ...
     .* exp(-s.^2 .* sin(a).^2) ...
   + sin(a)./s.^2 .* (1 + 2*s.^2 + e.*(1 - 2*s.^2 .* cos(2*a))) ...
     .* erf(s.*sin(a)) ...
   + (1-e)./s .* sqrt(pi) .* sin(a).^2 ...
     .* sqrt(Tw/T);



end