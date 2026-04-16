% Plot historgrams of mesh slopes
close all
clear
clc

% texture 1: rand 3
% texture 2: rand 1
% texture 3: rand 2

path = '/Users/antonburtsev_1/DocumentsLocal/DARPA/new_rand_specular/rand_1/';
name = 'anglesspec_70_';
cases = [0.1 0.5 1 2 3 5 7 9 10];
% cases = [0.1 0.5 1 2 3 5 7];


path = '/Users/antonburtsev_1/DocumentsLocal/DARPA/AFM_drag/Surf_slopes/';
name = 'Ge_Si_AFM_';
cases = [1 2];

% path = '/Users/antonburtsev/Documents/University/AUSTIN/DARPA/-Tools/rand_slopes_70_new/';
% name = 'rms';
% % cases = [0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.75 1];
% cases = [0.01 0.1 0.2 0.3 0.4 0.5 0.75 1];

% path = '.';
% name = 'zoom_';
% % cases = [10 20 5 2.5];
% cases = [2.5 5 10 20];

n = length(cases);

for i = 1:n
    
A = load([path,'/',name,num2str(cases(i)),'.txt']);
% A = load([path,'/',name,num2str(cases(i)),'nm.txt']);

theta_y = rad2deg(A(:,2)); % direction angle of normal from the +y axis
azim_y = rad2deg(A(:,4));  % azimuth angle from +y axis (projection of surface normal onto x-y plane)



maxval = max(theta_y);
nn = find(theta_y==maxval);
theta_y(nn) = [];

maxval = max(azim_y);
nn = find(azim_y==maxval);
azim_y(nn) = [];



meanslope(i) = mean(theta_y);
modeslope(i) = mode(theta_y);

figure(2)
% histogram(theta_y,'DisplayName',[name,' ',num2str(cases(i))]); hold on
h = histogram(theta_y,'DisplayStyle','stairs','LineWidth',2,'DisplayName',[name,' ',num2str(cases(i))]); hold on

nn = find(h.Values==max(h.Values));
binCenters = (h.BinEdges(1:end-1) + h.BinEdges(2:end)) / 2;
peak(i) =  binCenters(nn);

mean_angle(i) =  mean(theta_y);


figure(3)
% histogram(azim_y,'DisplayName',[name,' ',num2str(cases(i))]); hold on
h = histogram(azim_y,'DisplayStyle','stairs','LineWidth',2,'DisplayName',[name,' ',num2str(cases(i))],'Normalization', 'probability'); hold on

% nn = find(h.Values==max(h.Values));
binCenters = (h.BinEdges(1:end-1) + h.BinEdges(2:end)) / 2;
% peak(i) =  binCenters(nn);

[fwhmx(i)] = fwhm(h.Values,binCenters);

end
%%
figure(2)
xlabel({'$\theta_y ~(^\circ)$'},'Interpreter','latex','FontSize',20)
ylabel({'Count'},'Interpreter','latex','FontSize',20)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width4
set(gca,'TickLabelInterpreter','latex')

xlim([0 90])
xticks([0:10:90])

figure(3)
xlabel({'$\alpha_y ~(^\circ)$'},'Interpreter','latex','FontSize',20)
ylabel({'Count'},'Interpreter','latex','FontSize',20)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width4
set(gca,'TickLabelInterpreter','latex')

xlim([-90 90])
xticks([-100:20:-20,0,20:20:90])
set(gca,'xMinorTick','on')


% %% Full width at half maximum
% data = h.Values;
% x = binCenters;
% halfMax = (min(data) + max(data)) / 2;
% % Find where the data first drops below half the max.
% index1 = find(data >= halfMax, 1, 'first');
% % Find where the data last rises above half the max.
% index2 = find(data >= halfMax, 1, 'last');
% fwhm = index2-index1 + 1; % FWHM in indexes.
% % OR, if you have an x vector
% fwhmx = x(index2) - x(index1);
% 
% fwhmx
%%
function [fwhmx] = fwhm(data,x)
%UNTITLED Full width at half maximum
%   Detailed explanation goes here

halfMax = (min(data) + max(data)) / 2;
% Find where the data first drops below half the max.
index1 = find(data >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(data >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in indexes.
% OR, if you have an x vector
fwhmx = x(index2) - x(index1);

end
