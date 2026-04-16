close all
clear
clc

path = '/Users/antonburtsev_1/DocumentsLocal/DARPA/AFM/';
file = 'Ge_Si_ExpUCB_5x5.txt';
% file = 'Ge_Si_ExpUCB_100x100.txt';

% file = 'Ge_Si_100x100.txt';
% file = 'Si_100x100.txt';

% file = 'Au_Si_Exp_100x100.txt';
% file = 'Au_Si_Exp_5x5.txt';

% file = 'Si_0deg_5x5.txt';

% file = 'Ge_Si_UCB_5x5.txt';
% file = 'Ge_Si_ExpUCB_100x100.txt';
% file = 'Ge_Si_ExpUCB_2_100x100.txt';

% path = '/Users/antonburtsev_1/DocumentsLocal/DARPA/AFM/long_exposure/';
% file = 'Ge_Si_Exp_1h_5x5.txt';
% file = 'Ge_Si_5x5.txt';
% file = 'Kapton_5x5.txt';
% file = 'Kapton_Exp_1h_5x5.txt';
% file = 'Kapton_Exp_1h_100x100.txt';

% Plot region size
% xmin = 20; %20;
% xmax = 40; %90;
% ymin = 65; %0;
% ymax = 85; %65;

% xmin = 47.5;
% xmax = 52.5;
% ymin = 47.5;
% ymax = 52.5;

xmin = 15;
xmax = 500;
ymin = -500;
ymax = 500;

% xmin = 20;
% xmax = 90;
% ymin = 0;
% ymax = 100;

% xmin = 0.5;
% xmax = 1.5;
% ymin = 0.5;
% ymax = 1.5;

% xmin = 1;
% xmax = 2;
% ymin = 3;
% ymax = 4;

% xmin = 2.5;
% xmax = 3.5;
% ymin = 2;
% ymax = 3;

%%

[nx,ny] = get_size([path,file]);
fprintf('Data Size: %d x %d\n', nx, ny);

zres = 0.001; % AFM machine resolution 20A = 2nm
% z_lim = [-0.5 0.5];
z_lim = [-1 1];
% z_lim = [-10 10];

A = importdata([path,file],'\t',7);

x = A.data(:,1);
y = A.data(:,2);
z = A.data(:,3);

X = reshape(x,nx,ny);
Y = reshape(y,nx,ny);
Z = reshape(z,nx,ny);

global_zbar = mean(z); % global mean
m = and(and(x>xmin,x<xmax),and(y>ymin,x<ymax)); % logial index of plot range

yvec = Y(1,:);
xvec = X(:,1);
% Indices of plot range
mx = find(and(xvec>xmin,xvec<xmax));
my = find(and(yvec>ymin,yvec<ymax));

zbar = mean(z(m)); % mean in plot region


global_zrms = rms(z-global_zbar);       % gloabl rms
zrms = rms(z(m)-zbar);           % plot range rms

Zorig = Z;
% Remove tilt
Znew = remove_tilt(X,Y,Z);
Znew_local = remove_tilt(X(mx,my),Y(mx,my),Z(mx,my));

%% Plot Origianl Z as is
figure(1)
surf(X,Y,Zorig)
shading interp
view(0,90)
colormap(copper)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width
set(gca,'TickLabelInterpreter','latex')
xlabel({'$x (\mu m)$'},'Interpreter','latex','FontSize',20)
ylabel({'$y (\mu m)$'},'Interpreter','latex','FontSize',20)
c1 = colorbar('FontSize',16);
c1.Location = 'eastoutside';
c1.Label.String = '$z (nm)$';
c1.Label.Interpreter = 'latex';
c1.TickLabelInterpreter = 'latex';
c1.LineWidth = 1.5;
c1.Label.Rotation = 90;
c1.Label.FontSize = 26;
set(gca, 'Layer', 'top'); grid off

% data = z;
% zlim = [-max(zres,abs(min(data))) max(zres,abs(max(data)))];
% caxis(zlim)

%% Plot Z - global Zmean
figure(2)
% surf(X,Y,Z-global_zbar)
pcolor(X,Y,Znew)
shading interp
view(0,90)
colormap(copper)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width
set(gca,'TickLabelInterpreter','latex')
xlabel({'$x (\mu m)$'},'Interpreter','latex','FontSize',20)
ylabel({'$y (\mu m)$'},'Interpreter','latex','FontSize',20)
% title(['rms = ',num2str(global_zrms),' nm'],'FontSize',20)
title(['rms = ',sprintf('%.3f',global_zrms),' nm'],'FontSize',18)

c2 = colorbar('FontSize',16);
c2.Location = 'eastoutside';
c2.Label.String = '$z (nm)$';
c2.Label.Interpreter = 'latex';
c2.TickLabelInterpreter = 'latex';
c2.LineWidth = 1.5;
c2.Label.Rotation = 90;
c2.Label.FontSize = 26;
set(gca, 'Layer', 'top'); grid off

% data = z-zbar;
% zlim = [-max(zres,abs(min(data))) max(zres,abs(max(data)))];
caxis(z_lim)

r = rectangle('Position',[xvec(mx(1)) yvec(my(1)) xvec(mx(end))-xvec(mx(1)) yvec(my(end))-yvec(my(1))]);
r.LineWidth = 1.5;
r.LineStyle = '-';
r.EdgeColor = [1 1 1];
%% Plot Z in region - Zmean in that region
figure(3)
% surf(X(mx,my),Y(mx,my),Znew_local)
pcolor(X(mx,my),Y(mx,my),Znew_local)
shading interp
view(0,90)
colormap(copper)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width
% title(['rms = ',num2str(zrms),' nm'],'FontSize',18)
title(['rms = ',sprintf('%.3f',zrms),' nm'],'FontSize',18)
set(gca,'TickLabelInterpreter','latex')
xlabel({'$x (\mu m)$'},'Interpreter','latex','FontSize',20)
ylabel({'$y (\mu m)$'},'Interpreter','latex','FontSize',20)
c3 = colorbar('FontSize',18);
c3.Location = 'eastoutside';
c3.Label.String = '$z (nm)$';
c3.Label.Interpreter = 'latex';
c3.TickLabelInterpreter = 'latex';
c3.LineWidth = 1.5;
c3.Label.Rotation = 90;
c3.Label.FontSize = 26;
set(gca, 'Layer', 'top'); grid off

axis tight

% data = z-zbar;
% zlim = [-max(zres,abs(min(data))) max(zres,abs(max(data)))];
caxis(z_lim)

%%
% Exposed
% an = 1;
% xmin = 1.5;
% xmax = 2.5;
% ymin = 1;
% ymax = 2;
% an = 2;
% xmin = 1;
% xmax = 2;
% ymin = 3.5;
% ymax = 4.5;
% an = 3;
% xmin = 3.99;
% xmax = 4.99;
% ymin = 3;
% ymax = 4;


an = 1;
xmin = 1;
xmax = 2;
ymin = 1;
ymax = 2;
an = 2;
xmin = 2.5;
xmax = 4.5;
ymin = 2.5;
ymax = 4.5;
an = 3;
xmin = 2.8;
xmax = 3.8;
ymin = 0.2;
ymax = 1.2;

% Unexposed
% an = 1;
% xmin = 0.5;
% xmax = 1.5;
% ymin = 0.5;
% ymax = 1.5;

% an = 2;
% xmin = 1;
% xmax = 2;
% ymin = 3;
% ymax = 4;

% an = 3;
% xmin = 2.5;
% xmax = 3.5;
% ymin = 2;
% ymax = 3;

% an = 1;
% xmin = 0.1;
% xmax = 1.1;
% ymin = 2.6;
% ymax = 3.6;
% 
% an = 2;
% xmin = 3;
% xmax = 4;
% ymin = 3.5;
% ymax = 4.5;
% % % 
% an = 3;
% xmin = 2.5;
% xmax = 4.5;
% ymin = 0.5;
% ymax = 2.5;


mx = find(and(xvec>xmin,xvec<xmax));
my = find(and(yvec>ymin,yvec<ymax));
xw = xvec(mx(end))-xvec(mx(1));
yw = yvec(my(end))-yvec(my(1));
r = rectangle('Position',[xvec(mx(1)) yvec(my(1)) xw yw]);
r.LineWidth = 1.5;
r.LineStyle = '-';
r.EdgeColor = [1 1 1];

text(xvec(mx(1))+0.1*xw ,yvec(my(end))-0.2*yw,num2str(an),'FontSize',18,'Color','w')

xlim([0 5]);ylim([0 5]);c3.Ticks = [-1:0.2:1];
axis square
% set(gca,'xticks',[0:1:5])

% m = and(and(x>xmin,x<xmax),and(y>ymin,y<ymax)); % logial index of plot range
m = and(and(X>xmin,X<xmax),and(Y>ymin,Y<ymax)); % logial index of plot range
zrms = rms(Znew_local(m));    
zrms % box rms


% text(xvec(mx(1))+0.1*xw ,yvec(my(end))+0.2*yw,num2str(zrms),'FontSize',18,'Color','w')
%%
function [nx,ny] = get_size(path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(path, 'r'); % Open the file
line = fgetl(fid); % Read the first line
fclose(fid); % Close the file
% Extract numbers using regexp
nums = regexp(line, '\d+', 'match');
nums = str2double(nums);
% Assign to variables
nx = nums(1);
ny = nums(2);
end

function [Znew] = remove_tilt(X,Y,Z)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Remove linear tilt
[nx, ny] = size(Z);
A = [X(:), Y(:), ones(numel(Z),1)];
coeff = A \ Z(:);   % least-squares fit
Zfit = reshape(A * coeff, nx, ny);
Znew = Z - Zfit;
end

