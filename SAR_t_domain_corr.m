% ================== SAR Time-Domain Correlation (Simulated Targets) ==================
clearvars; close all; clc;

% ---------- 0) Basic Parameters (modify for tasks a–d) ----------
Fs   = 50e3;
Tp   = 5e-3;                 % Up-chirp duration
f1   = 2.4e9; f2 = 2.5e9;    % Frequency band (use 5.4–5.5 GHz for part d)
BW   = f2 - f1;
fc   = (f1 + f2)/2;
c    = 3e8;
N    = round(Fs*Tp);         % Samples per chirp (N = Fs*Tp)

lambda  = c/fc;
spacing = lambda/2;

% Rail positions (meters)
xp = -3.000 : spacing : 3.000;   % cross-range radar track
yp = 0;

% Target positions (meters)
x_pos = [-2.050  0.000  2.500];
y_pos = [0 0 0];
z_pos = [10 20 30];

% Time array (row vector)
t = linspace(0, Tp, N);         % 1 x N

% ---------- 1) Generate Simulated Data using function ----------
% Expect SAR_Data to be: numPositions x 1 x N (or numPositions x channels x N)
SAR_Data = SAR_Data_Generation(x_pos, y_pos, z_pos, xp, yp, N, t, fc, BW, Tp);

% ---------- 2) Define z-axis (down-range) ----------
dR   = c/(2*BW);                 % Range resolution ΔR
Rmax = N * dR / 2;               % Maximum unambiguous range
z_range = linspace(0, Rmax, N);  % 1 x N

% ---------- 3) Define x-axis (cross-range) ----------
x_range = xp;                    % 1 x M
M = numel(xp);

% ---------- 4) Initialize Image Matrix ----------
img = zeros(numel(x_range), numel(z_range)); % M x N

% ---------- 5) Prepare measurement matrix (M x N) ----------
% Adjust extraction depending on the exact SAR_Data shape:
% If SAR_Data is (M x 1 x N) then squeeze gives M x N:
meas = squeeze(SAR_Data(:,1,:));     % should give M x N
% If your SAR_Data is M x N directly, use meas = SAR_Data;

% Ensure meas orientation: rows = radar positions, cols = time samples
if size(meas,1) ~= M || size(meas,2) ~= N
    error('Unexpected meas size: expected [%d x %d], got [%d x %d].', M, N, size(meas,1), size(meas,2));
end

% Convert t to row vector (ensure 1 x N) and meas is M x N
t = reshape(t, 1, []);   % 1 x N

% ---------- 6) Time-Domain Correlation (pixel-wise focusing) ----------
tic;
for ix = 1:M
    x_img = x_range(ix);                  % scalar
    % cross-range distance vector from pixel to all radar pos (1 x M)
    dx = x_img - xp;                      % 1 x M
    % absolute not necessary inside sqrt, but keep sign-consistent:
    r_x = abs(dx);                        % 1 x M

    for iz = 1:N
        z_img = z_range(iz);              % scalar
        % slant range vector (1 x M)
        r = sqrt(r_x.^2 + z_img^2);      % 1 x M
        % two-way delays as column vector (M x 1)
        td = (2 * r / c).';              % make it M x 1 (transpose)

        % Build phase kernel (M x N) using implicit expansion:
        % phase(m,n) = exp(j*2*pi*(fc*td(m) + (BW/Tp)*td(m) * t(n)))
        phase = exp(1j * 2 * pi * ( fc * td + (BW / Tp) * (td .* t) )); % M x N

        % Correlation: sum over radar positions and time
        pixel = sum(sum(meas .* phase, 2));    % sum over rows then cols -> scalar
        img(ix, iz) = abs(pixel);
    end
end
t_focus = toc;

% ---------- 7) Range-squared compensation ----------
% z_range is 1 x N; we want to multiply each column (range bin) by z^2:
img = img .* (ones(M,1) * (z_range .^ 2));   % M x N .* M x N

% ---------- 8) Display SAR Image ----------
figure;
imagesc(z_range, x_range, img ./ max(img(:)));
axis xy;
xlabel('Down-range z (m)');
ylabel('Cross-range x (m)');
title(sprintf('SAR TDC | Fs=%.0f kHz | fc=%.2f GHz | Rail length=%.1f m', Fs/1e3, fc/1e9, (max(xp)-min(xp))));
colormap parula; colorbar;

fprintf('Processing time: %.3f s\n', t_focus);
fprintf('Range resolution ΔR = %.2f m, Max range Rmax = %.1f m (N=%d)\n', dR, Rmax, N);

% ---------- Cross Range Resolution ----------
L = max(xp) - min(xp);             % aperture length
R_ref = mean(z_pos);               % reference range (average target depth)
delta_x = lambda * R_ref / (2 * L);
fprintf('Cross Range Resolution = %.2f m\n', delta_x);
