close all
clear
clc

B = readmatrix('/Users/antonburtsev_1/Documents/University/AUSTIN/DARPA/XPS_data/xps_data_2.xlsx','Sheet','before_exposure');
A = readmatrix('/Users/antonburtsev_1/Documents/University/AUSTIN/DARPA/XPS_data/xps_data_2.xlsx','Sheet','after_exposure');

% species = { 'GeO_2',33.2;
%             'Ge_2O_3',32.4;
%             'GeO',31.4;
%             'GeO_2',30.3;
%             'Ge',29.6;
%     };
% species = { 'GeO_2',32.9;
%             'Ge_2O_3',32.3;
%             'GeO',31.5;
%             'GeO_2',29.3;
%             'Ge',28.7;
%     };
% species = { 'GeO_2',33.2;
%             'Ge_2O_3',32.4;
%             'GeO',31.4;
%             'GeO_2',30.3;
%             'Ge',28.7;
%     };
% n = size(species,1);

%% Before Exposure

species = { '32 eV',32;
            '29.8 eV',29.8;
            '28.9 eV',28.9;
            '28.3 eV',28.3;
    };
n = size(species,1);

norm_val = max(B(:,8));
B(:,2:end) = B(:,2:end)./norm_val;

figure(1)
t = tiledlayout(5,1,"TileSpacing","none");
nexttile([4 1]);
scatter(B(:,1),B(:,8),10,'ko','LineWidth',2,'DisplayName','Data'); hold on
plot(B(:,1),B(:,2),'r-','LineWidth',2,'DisplayName','Fit'); hold on
set(gca,'ColorOrderIndex',1)
plot(B(:,1),B(:,6),'LineWidth',2,'DisplayName','Ge^{0}'); hold on
plot(B(:,1),B(:,5),'LineWidth',2,'DisplayName','Ge^{1+}'); hold on
plot(B(:,1),B(:,3),'LineWidth',2,'DisplayName','Ge^{2+}'); hold on
plot(B(:,1),B(:,4),'LineWidth',2,'DisplayName','Ge^{3+}'); hold on
plot(B(:,1),B(:,7),'LineWidth',2,'DisplayName','Ge^{4+}'); hold on
set(gca,'xDir','reverse')
set(gca,'xTickLabel','')
set(gca,'xMinorTick','on')
set(gca,'yTick',[0.2:0.2:1])
xlim([27 35])
legend('Location','northwest','Box','off')


yl = ylim;
% ylim([yl(1) yl(2)+0.25*yl(2)])
ylim([yl(1) 1.3])
yl = ylim;

fprintf('After:\n')
pk = find(B(:,7)==max(B(:,7))); fprintf(['Ge^{4+} ',num2str(B(pk,1)),' eV\n'])
pk = find(B(:,4)==max(B(:,4))); fprintf(['Ge^{3+} ',num2str(B(pk,1)),' eV\n'])
pk = find(B(:,3)==max(B(:,3))); fprintf(['Ge^{2+} ',num2str(B(pk,1)),' eV\n'])
pk = find(B(:,5)==max(B(:,5))); fprintf(['Ge^{1+} ',num2str(B(pk,1)),' eV\n'])
pk = find(B(:,6)==max(B(:,6))); fprintf(['Ge^{0} ',num2str(B(pk,1)),' eV\n'])

for i=1:n
    % plot([cell2mat(species(i,2)) cell2mat(species(i,2))], [yl(1) yl(2)],'k--','LineWidth',1.5,'HandleVisibility','off')
    % text(cell2mat(species(i,2)), yl(2)-0.1*(1-mod(i,2))*yl(2)-0.05*yl(2), ['$\leftarrow\mathrm{',char(species(i,1)),'}$'],'FontSize',18,'Interpreter','latex')
    plot([cell2mat(species(i,2)) cell2mat(species(i,2))], [yl(1) yl(2)-0.1*(1-mod(i,2))*yl(2)-0.08*yl(2)],'k--','LineWidth',1.5,'HandleVisibility','off')
    text(cell2mat(species(i,2)), yl(2)-0.1*(1-mod(i,2))*yl(2)-0.05*yl(2), ['$\mathrm{',char(species(i,1)),'}$'],'FontSize',18,'Interpreter','latex','HorizontalAlignment','center')
end

ylabel({'Intensity  (Arb. units)'},'Interpreter','latex','FontSize',20)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width4
set(gca,'TickLabelInterpreter','latex')
xlim([27 35])

nexttile
plot(B(:,1),B(:,8)-B(:,2),'k-','LineWidth',2); hold on
set(gca,'xDir','reverse')
set(gca,'xMinorTick','on')
ylim([-0.015 0.015])

xlabel({'Binding energy (eV)'},'Interpreter','latex','FontSize',20)
ylabel({'Residual'},'Interpreter','latex','FontSize',20)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width4
set(gca,'TickLabelInterpreter','latex')
xlim([27 35])
%% After Exposure

species = { '32.4 eV',32.4;
            '31.8 eV',31.8;
            '30.7 eV',30.7;
            '28.7 eV',28.7;
            '28.1 eV',28.1;
    };
n = size(species,1);

norm_val = max(A(:,9));
A(:,2:end) = A(:,2:end)./norm_val;

figure(2)
t = tiledlayout(5,1,"TileSpacing","none");

nexttile([4 1]);
scatter(A(:,1),A(:,9),10,'ko','LineWidth',2,'DisplayName','Data'); hold on
plot(A(:,1),A(:,2),'r-','LineWidth',2,'DisplayName','Fit'); hold on
set(gca,'ColorOrderIndex',1)
plot(A(:,1),A(:,5),'LineWidth',2,'DisplayName','Ge^{0}'); hold on
plot(A(:,1),A(:,6),'LineWidth',2,'DisplayName','Ge^{1+}'); hold on
plot(A(:,1),A(:,3),'LineWidth',2,'DisplayName','Ge^{2+}'); hold on
plot(A(:,1),A(:,4),'LineWidth',2,'DisplayName','Ge^{3+}'); hold on
plot(A(:,1),A(:,7),'LineWidth',2,'DisplayName','Ge^{4+}'); hold on

set(gca,'xDir','reverse')
set(gca,'xTickLabel','')
set(gca,'xMinorTick','on')
set(gca,'yTick',[0.2:0.2:1])
xlim([27 35])
legend('Location','northwest','Box','off')

yl = ylim;
% ylim([yl(1) yl(2)+0.1*yl(2)])
ylim([yl(1) 1.3])
yl = ylim;

fprintf('After:\n')
pk = find(A(:,7)==max(A(:,7))); fprintf(['Ge^{4+} ',num2str(A(pk,1)),' eV\n'])
pk = find(A(:,4)==max(A(:,4))); fprintf(['Ge^{3+} ',num2str(A(pk,1)),' eV\n'])
pk = find(A(:,3)==max(A(:,3))); fprintf(['Ge^{2+} ',num2str(A(pk,1)),' eV\n'])
pk = find(A(:,6)==max(A(:,6))); fprintf(['Ge^{1+} ',num2str(A(pk,1)),' eV\n'])
pk = find(A(:,5)==max(A(:,5))); fprintf(['Ge^{0} ',num2str(A(pk,1)),' eV\n'])

for i=1:n
    % plot([cell2mat(species(i,2)) cell2mat(species(i,2))], [yl(1) yl(2)],'k--','LineWidth',1.5,'HandleVisibility','off')
    % text(cell2mat(species(i,2)), yl(2)-0.1*(1-mod(i,2))*yl(2)-0.05*yl(2), ['$\leftarrow\mathrm{',char(species(i,1)),'}$'],'FontSize',18,'Interpreter','latex')
    plot([cell2mat(species(i,2)) cell2mat(species(i,2))], [yl(1) yl(2)-0.1*(1-mod(i,2))*yl(2)-0.08*yl(2)],'k--','LineWidth',1.5,'HandleVisibility','off')
    text(cell2mat(species(i,2)), yl(2)-0.1*(1-mod(i,2))*yl(2)-0.05*yl(2), ['$\mathrm{',char(species(i,1)),'}$'],'FontSize',18,'Interpreter','latex','HorizontalAlignment','center')
end

ylabel({'Intensity (Arb. units)'},'Interpreter','latex','FontSize',20)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width4
set(gca,'TickLabelInterpreter','latex')
xlim([27 35])

nexttile
plot(A(:,1),A(:,9)-A(:,2),'k-','LineWidth',2); hold on
set(gca,'xDir','reverse')
set(gca,'xMinorTick','on')
ylim([-0.015 0.015])
% set(gca,'yTick',[-0.01,0,0.01])

xlabel({'Binding energy (eV)'},'Interpreter','latex','FontSize',20)
ylabel({'Residual'},'Interpreter','latex','FontSize',20)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width4
set(gca,'TickLabelInterpreter','latex')
xlim([27 35])