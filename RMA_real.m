clearvars -except sync backscatter Fs
clc
close all

% Read file
[data, Fs] = audioread('SAR_25sep_edited.wav');
if size(data,2) == 1
    error('Input audio must have two channels: backscatter and sync (square).');
end


sync = data(:, 2);
backscatter = data(:, 1);
clear data
%
% figure;
% 
% subplot(2,1,1)   % 2 rows, 1 column, first plot
% plot(backscatter);
% ylabel('Amplitude');
% xlabel('Data Sample Number');
% title('Radar backscatter data');
% 
% subplot(2,1,2)   % 2 rows, 1 column, second plot
% plot(sync);
% ylabel('Amplitude');
% xlabel('Data Sample Number');
% title('Sync data');




%%
% Make sync binary using a small adaptive threshold
thr_find_pos = 0.05;
thr = 0;
sync_data = sync > thr_find_pos;


% Expected number of samples per up-chirp
Trp = 0.3;              % range profile time in seconds (e.g. 275 ms)
Nrp = round(Trp * Fs);  % samples per range profile
guard = Nrp;
Tp = 0.005;

% Detect sync pulse edges
d = diff([0; sync_data; 0]);
starts = find(d == 1); % When each chirp starts
ends   = find(d == -1) - 1; % When each chirp ends

% Remove chirps that are too short or too long
for i = 1:length(starts)
    chirp_length = ends(i) - starts(i);
    if chirp_length < 0.5*Tp*Fs || chirp_length > 2*Tp*Fs
        starts(i) = 0; ends(i) = 0;
    end
end
starts(starts == 0) = []; ends(ends == 0) = [];

% Remove chirps that are too close together
for i = 1:length(starts)-1
    spacing_length = starts(i+1) - ends(i);
    if spacing_length < 0.5*Tp*Fs || spacing_length > 2*Tp*Fs
        starts(i) = 0; ends(i) = 0;
    end
end
starts(starts == 0) = []; ends(ends == 0) = [];


chirp_times = ends - starts;
silence_after = [starts(2:end); length(sync_data)]-ends(1:end); % How much silence is after each chirp
silence_before = [inf; silence_after(1:end-1)]; % How much silence is before each chirp
long_silences = silence_before > Nrp;

pos_starts = starts .* long_silences;
pos_starts(pos_starts == 0) = []; % Start time of each chirp that is preceeded by a long silence

sync_data = sync > thr;

% test123 = zeros(1, length(sync));
% test123(pos_starts) = 2;
% len_over_10 = floor(length(sync)/10);
% for i = 0:9
%     idx = len_over_10*i+1:len_over_10*(i+1);
%     figure
%     plot(sync(idx)); hold on; plot(test123(idx));
% end
% plot(sync); hold on; plot(test123);
% plot(sync(round(end*0.3):round(end*0.5)));hold on; plot(test123(round(end*0.3):round(end*0.5)))
% plot(sync(1:pos_starts(1)+Nrp*2));hold on; plot(test123(1:pos_starts(1)+Nrp*2))
% test321 = zeros(size(test123));
% test321(pos_starts(1)+Nrp) = 2;
% plot(test321(1:pos_starts(1)+Nrp*2))

% test123 = zeros(1, length(sync));
% for i = 1:length(pos_starts)
%     pos = pos_starts(i);
%     test123(pos + Nrp : pos + 2*Nrp-1) = sync(pos + Nrp : pos + 2*Nrp-1);
% end
% plot(sync(1:2000000))
% hold on
% plot(test123(1:2000000))
% title('Sync data')
% ylabel('Amplitude');
% xlabel('Data Sample Number');
% legend('Sync data', 'Parsed sync data')


%%
% Number of positions (one per sync pulse)
M = round(length(pos_starts));

% Keep only the first tablesworth of measurments
% third = floor(length(pos_starts)/3);
% pos_starts = pos_starts(1:third);
% M = length(pos_starts);

% Keep only an odd number of radar positions
if mod(M, 2) == 0
    M = M - 1;  % discard the last position if even
    pos_starts(end) = [];
end


%%
% Allocate matrices
DataMatrix = zeros(M, Nrp);
SyncMatrix = zeros(M, Nrp);

for m = 1:M
    s = pos_starts(m) + guard;  % start index after guard band
    if s+Nrp-1 <= length(backscatter)
        DataMatrix(m,:) = backscatter(s:s+Nrp-1);
        SyncMatrix(m,:) = sync_data(s:s+Nrp-1);
    else
        warning("Position %d truncated (not enough samples).", m);
    end
end

Tp = 0.005;
N = floor(Tp*Fs/2)*2;

% 1 = DC, N/2+1 = Nyquist

ProcessedUpChirps = zeros(M, N);
%ProcessedUpChirps = zeros(round(M/2), N);

for m = 1:M
    data_row = DataMatrix(m, :);
    sync_row = SyncMatrix(m, :);

    % Count number of up-chirps by detecting edges
    % d_sync = diff([0 sync_row 0]);
    % upchirp_starts = find(d_sync == 1);
    % upchirp_ends   = find(d_sync == -1) - 1;
    % numUpChirps = length(upchirp_starts);

    d_sync = diff(sync_row);
    upchirp_starts = find(d_sync == 1) + 1;
    upchirp_ends   = find(d_sync == -1);

    if upchirp_ends(1) < upchirp_starts(1)
        upchirp_ends(1) = [];
    end
    if upchirp_ends(end) < upchirp_starts(end)
        upchirp_starts(end) = [];
    end
    if length(upchirp_starts) ~= length(upchirp_ends)
        disp('FAULTY CHIRPS!')
        disp(upchirp_starts)
        disp(upchirp_ends)
    end
    numUpChirps = length(upchirp_starts);

    if numUpChirps == 0
        warning("No up-chirps found at position %d", m);
        continue
    end

    % Sum each integrated upchirp
    upchirp_sum = zeros(1, N);
    for i = 1:numUpChirps
        chirp_data = data_row(upchirp_starts(i):upchirp_ends(i));
        if length(chirp_data) < 10
            disp(length(chirp_data));
        end
        % Make the chirp data exactly N long
        t = linspace(0,1,length(chirp_data));
        tN = linspace(0,1,N);
        chirp_data = interp1(t, chirp_data, tN);

        upchirp_sum = upchirp_sum + chirp_data;
    end

    % Average over number of up-chirps
    avg_upchirp = upchirp_sum / numUpChirps;

    % Apply Hann window
    w = hann(length(avg_upchirp))';
    avg_upchirp = avg_upchirp .* w;

    % Hilbert transform
    avg_upchirp = hilbert(avg_upchirp);
    avg_upchirp(isnan(avg_upchirp)) = 1e-30;

    % Put the values into the matrix
    ProcessedUpChirps(m, :) = avg_upchirp;
%     if m == 5   % <-- choose position you want to visualize
%     figure;
% 
%     % Plot raw data with sync overlay
%     subplot(2,1,1);
%     plot(data_row); hold on;
%     yyaxis right
%     plot(sync_row, 'r'); % sync signal
%     yyaxis left
%     xlim([0 5000])
% 
%     % Mark chirp start/ends
%     for i = 1:numUpChirps
%         xline(upchirp_starts(i), 'g--');
%         xline(upchirp_ends(i), 'r--');
%     end
%     title(sprintf('Data Row with Sync (Position %d)', m));
%     xlabel('Sample Index'); ylabel('Amplitude');
%     grid on;
% 
%     % Plot averaged chirp
%     subplot(2,1,2);
%     plot(real(avg_upchirp), 'b', 'LineWidth', 2);
%     title(sprintf('Averaged Up-Chirp (Position %d)', m));
%     xlabel('Sample Index'); ylabel('Amplitude');
%     grid on;
% end

end


%%
% figure
% imagesc(ProcessedUpChirps);
% xlabel('Up-chirp fast time (samples)'); ylabel('Rail position'); 
% title('Integrated up-chirp data (rows=positions, cols=fast-time)');

%%
% Radar parameters
Tp = 5e-3;
f_start = 2.4e9; f_stop = 2.5e9;

% Fast time axis
t = linspace(0, Tp, N);

% Slow time / radar rail parameters
spacing = 0.06;           % step size (m)
L = spacing * size(ProcessedUpChirps,1);   % total rail length (m)

% Slow time axis / radar positions
xp = linspace(-L/2, L/2, size(ProcessedUpChirps,1));

% Constants
c = 299792458;              % Speed of light (m/s)
fc = (f_stop + f_start)/2;  % Carrier frequency (Hz), midpoint of 2.4-2.5 GHz
BW  = f_stop - f_start; % Bandwidth (Hz)

% Generate kt
kt = 2*pi*fc/c + (2*pi*BW/Tp)*(t/c);

% figure
% imagesc(kt, xp, angle(ProcessedUpChirps));
% xlabel('K_t (rad/m)');
% ylabel('Synthetic aperture position x_p (m)');
% title('Phase before along track FFT');

% Desired zero-padding size
zpad = 2049;

% Calculate number of zeros to add on top and bottom
pad_top    = floor((zpad - M)/2);
pad_bottom = ceil((zpad - M)/2);

% Add zero-padding symmetrically
ProcessedUpChirps_zp = [zeros(pad_top, N); ProcessedUpChirps; zeros(pad_bottom, N)];

% Current size
M = size(ProcessedUpChirps_zp, 1); % M = 2049

% Updated L and xp
L = spacing * M  ;   % total rail length (m)
xp = linspace(-L/2, L/2, M);

% Generate kx_2D
delta_kx = pi/max(xp); % 2*pi/(2*max(xp))???
kx = delta_kx * ( -(size(xp, 2)-1)/2 : (size(xp, 2)-1)/2 );
kx_2D = repmat(kx.', [1,N]);

% Generate kt_2D
kt_2D = repmat(kt, [size(xp,2),1]);

%Generate kz_2D
kz_2D = sqrt(4*(kt_2D.^2) - kx_2D.^2);

% Generate uniform kz
kz_1D_uni = linspace(min(min(kz_2D)), max(max(kz_2D)), size(kz_2D,2));
kz_2D_uni = repmat(kz_1D_uni,[size(xp,2),1]);

% figure
% scatter(kx_2D(1:10:end), kz_2D(1:10:end));
% hold on
% scatter(kx_2D(1:10:end), kz_2D_uni(1:10:end));
% title('K_z comaprison');

% Perform 1-D FFT
S_B = fftshift(fft(ProcessedUpChirps_zp,[],1),1);

% figure
% imagesc(kt, kx, 20*log10(abs(S_B)), [-18, 20.5]);
% xlabel('K_t (rad/m)')
% ylabel('K_x (rad/m)')
% title('Magnitude after along tarck FFT');
% colorbar;
% 
% figure
% imagesc(kt, kx, angle(S_B));
% xlabel('K_t (rad/m)')
% ylabel('K_x (rad/m)')
% title('Phase after along tarck FFT');
% colorbar;
%%

% Stolt interpolation
S_B2 = zeros(size(xp, 2), N);
for i = 1:size(xp, 2)
    S_B2(i,:) = interp1(kz_2D(i,:), S_B(i,:), kz_1D_uni);
end
S_B2(isnan(S_B2)) = 1e-30;

% figure
% imagesc(kz_1D_uni, kx, angle(S_B2))
% xlabel('K_z uniform (rad/m)')
% ylabel('K_x (rad/m)')
% title('Phase after along tarck FFT');
% colorbar;

% 2-D IFFT
ft_1 = ifft(S_B2,[],1);
ft_2 = ifft(ft_1,[],2);

% Create indices
delta_fz = c * (max(kz_1D_uni) - min(kz_1D_uni)) / (4*pi);
delta_kz = 4*pi*delta_fz/c;
Rmax = N*c/(2*delta_fz);
delta_x = 0.06; % ???
Rail_Rmax = size(ft_2, 1) * delta_x;

d_range_1 = 1; d_range_2 = floor(Rmax/4);
c_range_1 = -10; c_range_2 = 10;

d_index_1 = ceil((size(ft_2, 2)/Rmax)*d_range_1);
d_index_2 = ceil((size(ft_2, 2)/Rmax)*d_range_2);

c_index_1 = ceil((size(ft_2, 1)/Rail_Rmax)*(c_range_1+(Rail_Rmax/2)));
c_index_2 = ceil((size(ft_2, 1)/Rail_Rmax)*(c_range_2+(Rail_Rmax/2)));

% Flip and rotate
ft_2 = rot90(ft_2, 1);
ft_2 = fliplr(ft_2);

% Truncate
truncated_data_matrix = ft_2(d_index_1:d_index_2, c_index_1:c_index_2);

% Create downrange and crossrange
downrange = linspace(-1*d_range_1, - 1*d_range_2, size(truncated_data_matrix, 1));
crossrange = linspace(c_range_1, c_range_2, size(truncated_data_matrix, 2));

% Combine truncated data with downrange vector
truncated_data_matrix = (downrange(:).^2) .* truncated_data_matrix;

% Create image
image = 20*log10(abs(truncated_data_matrix));

% Plot image
figure
imagesc(crossrange, downrange, image, ...
    [max(max(image))-25, max(max(image))+0]);
colorbar;
xlabel('Crossrange (meter)');
ylabel('Downrange (meter)');
title('Final image');