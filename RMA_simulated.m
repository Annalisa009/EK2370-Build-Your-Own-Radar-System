close all

% Radar parameters
Fs = 50e3;
Tp = 5e-3;
f_start = 2.4e9; f_stop = 2.5e9;


% Data generation inputs
c = 299792458;
x_pos = [-2050e-3 0e-3 2500e-3]; y_pos = [0 0 0]; z_pos = [10 20 30];
fc = (f_start + f_stop) / 2;
lambda = c/fc; spacing = lambda/2;
xp = -3000e-3:spacing:3000e-3; yp = 0;
N = floor(Tp*Fs);
t = linspace(0, Tp, N);
BW = f_stop - f_start;

% Generate data
SAR_data = SAR_Data_Generation(x_pos, y_pos, z_pos, xp, yp, N, t, fc, BW, Tp);
data=squeeze(SAR_data);

% Generate kx_2D
delta_kx = pi/max(xp); % 2*pi/(2*max(xp))???
kx = delta_kx * ( -(size(xp, 2)-1)/2 : (size(xp, 2)-1)/2 );
kx_2D = repmat(kx.', [1,N]);

% Generate kt_2D
kt = 2*pi*fc/c + (2*pi*BW/Tp)*(t/c);
kt_2D = repmat(kt, [size(xp,2),1]);

% Generate kz_2D
kz_2D = sqrt(4*(kt_2D.^2) - kx_2D.^2);

% Generate uniform kz
kz_1D_uni = linspace(min(min(kz_2D)), max(max(kz_2D)), size(kz_2D,2));
kz_2D_uni = repmat(kz_1D_uni,[size(xp,2),1]);

% Plot the K-space both the nonuniformly and uniformly sampled kz for comparison
scatter(kx_2D(1:10:end), kz_2D(1:10:end)); hold on
scatter(kx_2D(1:10:end), kz_2D_uni(1:10:end));
xlabel('k_x');
ylabel('k_z');

% Cross range 1-D FFT with no zero padding
S_B = fftshift(fft(SAR_data,[],1),1);

% Stolt interpolation
S_B2 = zeros(size(xp, 2), N);
for i = 1:size(xp, 2)
    S_B2(i,:) = interp1(kz_2D(i,:), S_B(i,:), kz_1D_uni);
end
S_B2(isnan(S_B2)) = 1e-30;

% 2-D IFFT
ft_1 = ifft(S_B2,[],1);
ft_2 = ifft(ft_1,[],2);

% Generate the axis for the final image
delta_x = pi/max(kx); % 2*pi/(2*max(kx))???
x_axis = delta_x * ( -(size(xp, 2)-1)/2 : (size(xp, 2)-1)/2 );
delta_fz = c * (max(kz_1D_uni) - min(kz_1D_uni)) / (4*pi);
Rmax = N*c / (4*delta_fz);
z_axis = linspace(0, Rmax, N/2);

% Multiply each column, element by element by square of the range axis.
% Use lower half of the data
f_corrected = ft_2(:,1:size(ft_2,2)/2).* z_axis.^2;

figure
imagesc(z_axis, x_axis,(abs(f_corrected/max(f_corrected(:)))));
colorbar;
title('SAR Image');
xlabel('Z (m)');
ylabel('X (m)');


function SAR_Data_Time_Domain = SAR_Data_Generation(x_pos, y_pos, z_pos, xp, yp, N, t, fc, BW, Tp)

    c = 299792458; %(m/s) speed of light
    f = 1*ones(1, length(x_pos)); % reflectivity function for the targets
    n_targets = length(x_pos); % Total number of targets
    Tot_x = length(xp); % Total number of x positions/measurements
    Tot_y = length(yp); % Total number of y positions/measurements
    
    SAR_Data_Time_Domain = zeros(Tot_x,Tot_y,N);
    
    for nt = 1:n_targets
        for m = 1:length(xp)
            r = sqrt((x_pos(nt)-xp(m)).^2 + y_pos(nt).^2 + z_pos(nt).^2);
            tau = 2 * r / c;
            phase = cos((2*2*pi*fc*r/c)+(2*2*pi*BW*r*t/(c*Tp)));  % 1×N
            SAR_Data_Time_Domain(m, 1, :) = SAR_Data_Time_Domain(m, 1, :) + reshape(f(nt)*(1./(r.^2)).*phase, 1, 1, []);
        end
    end
end