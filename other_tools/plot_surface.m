close all
clear
clc

% % path = '/Users/antonburtsev/DocumentsLocal/DARPA/Profilometry/Exposed/';
% path = '/Users/antonburtsev/DocumentsLocal/DARPA/Profilometry/04_09_2025/';
% 
% % filename = '2.5x_Height';
% % filename = '5x_Height';
% % filename = '10x_Height';
% % filename = '20x_Height';
% % filename = '50x_Height';
% % filename = '150x_Height';
% filename = 'A_exposed_10x_Stitching_repeated_2_Height';


path = '/Volumes/T7/AUSTIN/DARPA/Profilometry/08_09_2025/Ge_Si_150x/';
filename = 'Height_Processed';

% path = '/Users/antonburtsev/DocumentsLocal/DARPA/Profilometry/08_09_2025/Si_150x/';
% filename = 'Height_processed';
% filename = 'single_Height_Processed';

% path = '/Users/antonburtsev/DocumentsLocal/DARPA/Profilometry/08_09_2025/Ge_Kapton_150x/';
% filename = 'Height_processed';

% path = '/Users/antonburtsev/DocumentsLocal/DARPA/Profilometry/08_09_2025/Ge_Kapton_Exposed_Minton_150x/';
% filename = 'Height_Processed';

% path = '/Users/antonburtsev/DocumentsLocal/DARPA/Profilometry/08_09_2025/Ge_Si_Exposed_Minton_150x/';
% filename = 'Height_Processed';

% path = '/Users/antonburtsev/DocumentsLocal/DARPA/Profilometry/Images_for_Fiaz/Before_exposure/';
% filename = 'Height_corrected';

process = 0;
kx = 10;  % number of low modes in x-direction
ky = 10;  % number of low modes in y-direction

%%
headerlines = 15;
[xyCalibration, header] = read_header([path,filename,'.csv'],headerlines);
data = readmatrix([path,filename,'.csv'], 'NumHeaderLines', headerlines);

% scale = 1e-6; % cover micrometers to meters
scale = 1;

delta = xyCalibration*scale;
nx=size(data,2);
ny=size(data,1);

lx = nx*delta;
ly = ny*delta;
x = linspace(-0.5*lx,0.5*lx,nx);
y = linspace(-0.5*ly,0.5*ly,ny);

min(x)
min(y)

[X,Y] = meshgrid(x,y);

data = data*scale;

[nr, nc] = find(isnan(data));
data(nr,nc) = 0;
X(nr,nc) = 0;
Y(nr,nc) = 0;


%%
Z = data;  % Your surface data (MxN matrix)


% Z = data - mean(mean(data));  % Your surface data (MxN matrix)
if process == 1
    % Step 1: Apply 2D DCT
    Z_dct = dct2(Z);
    
    % figure
    % imagesc(Z_dct)
    
    % Step 2: Zero out high-frequency components (keep only low-wavelength modes)
    % Let's say you want to keep only the top-left kx × ky block
    
    % Create a mask that keeps only low-frequency components
    [M, N] = size(Z_dct);
    mask = zeros(M, N);
    mask(1:kx, 1:ky) = 1;
    
    % Extract the low-frequency part
    Z_low = idct2(Z_dct .* mask);
    
    % Step 3: Subtract the low modes from the original to get the residual (roughness)
    Z_high = Z - Z_low;
    
    
    % Optional: Visualize
    figure;
    subplot(1,3,1); surf(X,Y,Z,'FaceColor',[.5 .5 .5],'FaceLighting','gouraud'); title('Original Surface'); shading interp; colorbar;
    view(0,90); axis tight
    % view(35,35)
    % zlim(1.0e-03*[0.1500    0.4000])
    subplot(1,3,2); surf(X,Y,Z_low,'FaceColor',[.5 .5 .5],'FaceLighting','gouraud'); title('Low-Wavelength Modes'); shading interp; colorbar;
    view(0,90); axis tight
    % view(35,35)
    % zlim(1.0e-03*[0.1500    0.4000])
    subplot(1,3,3); surf(X,Y,Z_high,'FaceColor',[.5 .5 .5],'FaceLighting','gouraud'); title('High-Frequency Residual'); shading interp; colorbar;
    view(0,90); axis tight
    % view(35,35)
    % zlim(1.0e-03*[0.1500    0.4000])
    
    
    % roughness in micrometers
    height_rms = rms((data-mean(mean(data)))./scale,"all")
    height_rms = rms(Z_high./scale,"all")
    

else
    Z_high = Z;
end
%% RMS slope
Zmean = mean(mean(Z_high));
Sa = mean(mean(Z_high));
Sq = sqrt(1/(nx*ny) * sum((Z_high-Zmean).^2,'all'));

[dZdx,dZdy] = gradient(Z_high*1e3,delta); % conver to nm
Area = (abs(min(x))+max(x))*1e3*(abs(min(y))+max(y))*1e3; % area in nm^2
% Sdq = sqrt(1/Area * (sum(dZdx.^2,'all')+sum(dZdx.^2,'all')));
Sdq = sqrt(1/(nx*ny) * (sum(dZdx.^2,'all')+sum(dZdy.^2,'all')));

spec = [Sa, Sq, Sdq];

%% Wavelength distribution
[lambda_x, PSD_x, lambda_y, PSD_y] = decompose_wavelengths_matlab(Z_high, delta, delta);

figure(1)
loglog(lambda_x,PSD_x,'LineWidth',1.5); hold on
loglog(lambda_y,PSD_y,'LineWidth',1.5); hold on

xlabel({'Wavelength ~$(\mu m)$'},'Interpreter','latex','FontSize',20)
ylabel({'PSD'},'Interpreter','latex','FontSize',20)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width4
set(gca,'TickLabelInterpreter','latex')

%% Plot Height
figure
surf(X,Y,Z_high*1e3,'FaceColor',[.5 .5 .5],'FaceLighting','gouraud'); shading interp; colorbar;
view(0,90); axis tight

xlabel({'$x~ (\mu m)$'},'Interpreter','latex','FontSize',20)
ylabel({'$z~ (\mu m)$'},'Interpreter','latex','FontSize',20)
zlabel({'$z~ (n m)$'},'Interpreter','latex','FontSize',20)
box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width4
set(gca,'TickLabelInterpreter','latex')
set(gca, 'Layer', 'top'); grid off
c = colorbar('FontSize',18);
c.Location = 'eastoutside';
c.Label.String = 'Height ($n m$)';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.LineWidth = 1.5;
c.Label.Rotation = 90;
c.Label.FontSize = 20;

%     caxis([-50 50])
%     c.Ticks = [-50:10:50];

% view(53,29)
% pbaspect([1 1 0.1])

colormap(bluewhitered)

%%
% writecell(header',[path,filename,'_corrected.csv'],'Delimiter',',')
% writematrix(Z_high./scale,[path,filename,'_corrected.csv'],'WriteMode','append')
%% RMS
R = rms_overlapping_boxes(Z_high, 20);
% R = rms_overlapping_boxes(data(570:570+90,512:512+90), 20);

% R = R(1390:2232,1855:2908);
% X = X(1390:2232,1855:2908);
% Y = Y(1390:2232,1855:2908);

figure
surf(X,Y,R*1e3,'FaceColor',[.5 .5 .5],'FaceLighting','gouraud')
% surf(X(570:570+90,512:512+90),Y(570:570+90,512:512+90),R,'FaceColor',[.5 .5 .5],'FaceLighting','gouraud')
colorbar;
shading interp
view(0,90)
axis equal
axis tight
set(gca, 'Layer', 'top');

caxis([0 0.1])
colormap(hot)

box on;
set(gca,'FontSize',18)      % Axis and labels fontsise
set(gca,'linewidth',1.5)    % Axis line width
set(gca,'TickLabelInterpreter','latex')
xlabel({'$x ~(\mu m)$'},'Interpreter','latex','FontSize',20)
ylabel({'$y ~(\mu m)$'},'Interpreter','latex','FontSize',20)
c = colorbar('FontSize',18);
c.Location = 'eastoutside';
% c.Label.String = 'RMS ($\mu m$)';
c.Label.String = 'RMS ($n m$)';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.LineWidth = 1.5;
c.Label.Rotation = 90;
c.Label.FontSize = 20;

ind = isnan(R);
mean_rms = mean(R(not(ind)))*1e3

%% RMS wavelength spectra

% [lambda, rms_per_lambda] = roughness_vs_wavelength(Z_high, lx);
% figure(6);
% loglog(lambda*1e6, rms_per_lambda*1e6, 'LineWidth', 2); hold on
% xlabel('\lambda [\mum]');
% ylabel('RMS Roughness per \lambda [\mum]');
% grid on;
% title('Surface Roughness Spectrum');

% Measure roughness along x and y
[lambda_x, rms_x, lambda_y, rms_y] = roughness_vs_wavelength_dir(Z_high, lx, ly);

% Plot results
figure(6);
loglog(lambda_x*1e6, rms_x(1:end-1)*1e9, 'b', 'LineWidth', 2); hold on;
loglog(lambda_y*1e6, rms_y(1:end-1)*1e9, 'r', 'LineWidth', 2);
xlabel('\lambda [\mum]');
ylabel('RMS Roughness per \lambda [nm]');
title('Directional Roughness Spectrum');
legend('Horizontal (x)', 'Vertical (y)');
grid on;


%%
function [xyCalibration, header] = read_header(filename,headerlines)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filename, 'r');

xyCalibration = NaN;  % Default if not found
maxLines = headerlines;

for i = 1:maxLines
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    
    header{i} = line;
    
    % Debug: show line number and content
    fprintf('Line %d: %s\n', i, line);

    if contains(line, 'XY calibration', 'IgnoreCase', true)
        parts = strsplit(line, '\t');  % Try tab first
        if numel(parts) < 2
            parts = strsplit(line);    % Try space split if tab fails
        end

        % Remove units like 'um' if present
        valueStr = regexprep(parts{2}, '[^\d\.\-eE+]', '');
        xyCalibration = str2double(valueStr);
%         break;
    end
end

fclose(fid);


end
%%
function R = rms_overlapping_boxes(H, w)
    % RMS_OVERLAPPING_BOXES computes the RMS of height values in overlapping w×w boxes
    %
    % Inputs:
    %   H - 2D matrix of height values
    %   w - box width (integer)
    %
    % Output:
    %   R - 2D matrix of RMS values (same size as H, NaN where window doesn't fit)

    [nRows, nCols] = size(H);
    R = nan(nRows, nCols);  % initialize output with NaNs

    % Precompute half window size
    half_w = floor(w / 2);

    % Loop over valid centers
    for i = 1+half_w : nRows-half_w
        for j = 1+half_w : nCols-half_w
            % Extract the local box
            box = H(i-half_w:i+half_w, j-half_w:j+half_w);
            
%             box = box - mean(box);

            % Compute RMS
            R(i,j) = sqrt(mean(box(:).^2));
        end
    end
end
%%
function [lambda_x, rms_x, lambda_y, rms_y] = roughness_vs_wavelength_dir(h, Lx, Ly)
    % Computes roughness vs wavelength along X and Y directions separately.
    %
    % INPUTS:
    %   h   -> height map [Ny x Nx] in meters
    %   Lx  -> physical size along x [m]
    %   Ly  -> physical size along y [m]
    %
    % OUTPUTS:
    %   lambda_x -> wavelengths along x [m]
    %   rms_x    -> RMS roughness contribution per λ along x [m]
    %   lambda_y -> wavelengths along y [m]
    %   rms_y    -> RMS roughness contribution per λ along y [m]
    
    [Ny, Nx] = size(h);

    % Remove mean height
    h = h - mean(h(:));

    %% ----- 1. FFT along X (horizontal spectrum) -----
    Hx = fft(h, Nx, 2); % FFT along rows (x-direction)
    Pxx = mean(abs(Hx).^2, 1) / Ny; % Average PSD along y

    % Wavenumbers along x
    kx = (2*pi/Lx) * [0:(Nx/2-1), -Nx/2:-1];
    kx = fftshift(kx);
    Pxx = fftshift(Pxx); % Align zero frequency in the center

    % Take only positive frequencies
    kx_pos = kx(Nx/2+1:end);
    Pxx_pos = Pxx(Nx/2+1:end);

    % Convert to wavelength
    lambda_x = 2*pi ./ kx_pos;

    % Convert PSD to RMS per wavelength bin
    dkx = abs(kx(2) - kx(1));
    rms_x = sqrt(Pxx_pos * dkx);

    %% ----- 2. FFT along Y (vertical spectrum) -----
    Hy = fft(h, Ny, 1); % FFT along columns (y-direction)
    Pyy = mean(abs(Hy).^2, 2) / Nx; % Average PSD along x
    Pyy = Pyy(:)'; % Row vector

    % Wavenumbers along y
    ky = (2*pi/Ly) * [0:(Ny/2-1), -Ny/2:-1];
    ky = fftshift(ky);
    Pyy = fftshift(Pyy);

    % Take only positive frequencies
    ky_pos = ky(Ny/2+1:end);
    Pyy_pos = Pyy(Ny/2+1:end);

    % Convert to wavelength
    lambda_y = 2*pi ./ ky_pos;

    % Convert PSD to RMS per wavelength bin
    dky = abs(ky(2) - ky(1));
    rms_y = sqrt(Pyy_pos * dky);

    %% ----- 3. Remove NaN or Inf -----
%     valid_x = ~isnan(rms_x) & ~isinf(rms_x);
%     lambda_x = lambda_x(valid_x);
%     rms_x = rms_x(valid_x);
% 
%     valid_y = ~isnan(rms_y) & ~isinf(rms_y);
%     lambda_y = lambda_y(valid_y);
%     rms_y = rms_y(valid_y);
end

%%
function [lambda, rms_per_lambda] = roughness_vs_wavelength(h, L)
    % Computes roughness vs wavelength from a height map using PSD analysis
    %
    % INPUTS:
    %   h  -> NxN matrix of height values [m]
    %   L  -> physical size of domain [m]
    %
    % OUTPUTS:
    %   lambda         -> vector of wavelengths [m]
    %   rms_per_lambda -> RMS roughness per wavelength bin [m]
    
    % Get grid size
    [Ny, Nx] = size(h);
    N = min(Nx, Ny);
    dx = L / N;
    
    % Remove mean height
    h = h - mean(h(:));
    
    % 2D FFT of height map
    H = fftshift(fft2(h));
    P = abs(H).^2 / (N^2);  % power spectrum
    
    % Frequency grid
    kx = (2*pi/L) * [0:(N/2-1), -N/2:-1];
    ky = kx;
    [KX, KY] = meshgrid(kx, ky);
    K = sqrt(KX.^2 + KY.^2);
    
    % Bin data radially in k-space
    kmax = max(K(:));
    nbins = floor(N/2);
    kbins = linspace(0, kmax, nbins+1);
    Pk = zeros(1, nbins);
    counts = zeros(1, nbins);
    
    for i = 1:nbins
        mask = (K >= kbins(i)) & (K < kbins(i+1));
        Pk(i) = sum(P(mask));
        counts(i) = sum(mask(:));
    end
    
    % Average PSD per ring
    Pk = Pk ./ counts;
    
    % Convert to wavelength
    kcenter = 0.5 * (kbins(1:end-1) + kbins(2:end));
    lambda = 2*pi ./ kcenter;
    
    % Convert PSD to RMS roughness per wavelength bin
    dk = diff(kbins);
    rms_per_lambda = sqrt(Pk .* dk);
    
    % Remove NaN or Inf (at zero frequency)
    valid = ~isnan(rms_per_lambda) & ~isinf(rms_per_lambda);
    lambda = lambda(valid);
    rms_per_lambda = rms_per_lambda(valid);
end
%%
% H : M-by-N matrix of heights
% dx : spacing in x (distance between columns)
% dy : spacing in y (distance between rows)
function [lambda_x, PSD_x, lambda_y, PSD_y] = decompose_wavelengths_matlab(H, dx, dy)
    [M, N] = size(H);

    % remove mean (and optionally linear trend)
    H0 = H - mean(H(:));

    % window sizes (1D window applied along transform dimension)
    win_x = hann(N)';   % row-wise (length N)
    win_y = hann(M);    % column-wise (length M)

    % --- Horizontal decomposition (along x): FFT of each row, average PSD across rows ---
    % apply window along each row
    Hrows = (H0 .* repmat(win_x, M, 1));
    Nfft_x = 2^nextpow2(N);
    Fx = fft(Hrows, Nfft_x, 2);        % FFT along dim 2 (columns)
    % one-sided indices
    kx = 0:floor(Nfft_x/2);
    Fx_one = Fx(:, kx+1);
    PSD_rows = (abs(Fx_one).^2) / (Nfft_x^2);      % power (not exact units; consistent)
    % average PSD across rows
    PSD_x = mean(PSD_rows, 1);
    % spatial frequency (cycles per unit length)
    fx = kx ./ (Nfft_x * dx);
    lambda_x = NaN(size(fx));
    nonzero = fx > 0;
    lambda_x(nonzero) = 1 ./ fx(nonzero);
    lambda_x(~nonzero) = Inf;  % DC -> infinite wavelength

    % --- Vertical decomposition (along y): FFT of each column, average PSD across columns ---
    Hcols = (H0 .* repmat(win_y, 1, N));
    Nfft_y = 2^nextpow2(M);
    Fy = fft(Hcols, Nfft_y, 1);        % FFT along dim 1 (rows)
    ky = 0:floor(Nfft_y/2);
    Fy_one = Fy(ky+1, :);
    PSD_cols = (abs(Fy_one).^2) / (Nfft_y^2);
    PSD_y = mean(PSD_cols, 2).';      % row vector
    fy = ky ./ (Nfft_y * dy);
    lambda_y = NaN(size(fy));
    nonzero = fy > 0;
    lambda_y(nonzero) = 1 ./ fy(nonzero);
    lambda_y(~nonzero) = Inf;

    % return as column vectors for convenience
    lambda_x = lambda_x(:);
    PSD_x = PSD_x(:);
    lambda_y = lambda_y(:);
    PSD_y = PSD_y(:);
end



