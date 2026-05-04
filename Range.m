clear all; close all; clc
% Read the file
[data, Fs] = audioread("fmcw_single.wav");

% Make time axis for plotting the waveforms
t = (0:size(data,1)-1)/Fs;

% Separate the data and Visualise
sync_data = data(:, 2);
back_scatter_data = data(:, 1);
%plot(t, sync_data)
%figure
%plot(t, back_scatter_data)
%xlabel('Time (s)'), ylabel('Amplitude'), title('Mono Waveform')

% Split the Sync data
threshold = 0;
sync_data = sync_data > threshold;

% define up chirp duration and samples
Tp = 0.005;
chirp_duration = round(Tp*Fs);

% find the data where sync pulse is one
combined_data = sync_data.*back_scatter_data;

% find the segments in the data where sync pulse is one and concatenate them into a matrix
d = diff([0 transpose(sync_data) 0]);
starts = find(d == 1);
ends   = find(d == -1) - 1;
num_segments = length(starts);
seg_lengths = ends - starts + 1;
max_length = max(seg_lengths);
segments = zeros(num_segments, max_length);

for k = 1:num_segments
    seg = combined_data(starts(k): min(ends(k), starts(k) + chirp_duration+1));
    segments(k, 1:length(seg)) = seg;
end

% MS clutter rejection
col_mean = mean(segments, 1);
segments = segments - col_mean;

% generating the time axis
time = starts/Fs;

% performing MTI
data_2MTI = segments(2:size(segments,1),:) - segments(1:size(segments,1)-1,:);
data_3MTI = segments(3:size(segments,1),:) - 2*segments(2:size(segments,1)-1,:) ...
    + segments(1:size(segments,1)-2,:);

% Performing ifft
N = chirp_duration;
Y1= abs(ifft(segments, 4*N, 2));
Y2= abs(ifft(data_2MTI, 4*N, 2));
Y3= abs(ifft(data_3MTI, 4*N, 2));

% converting to db scale
Y1_dB = 20 * log10(Y1);
Y2_dB = 20 * log10(Y2);
Y3_dB = 20 * log10(Y3);

% taking the lower half of the matrix
half_size = floor(size(Y1_dB, 2)/2);
Y1_dB_lower = Y1_dB(:, 1:half_size);
Y2_dB_lower = Y2_dB(:, 1:half_size);
Y3_dB_lower = Y3_dB(:, 1:half_size);

% Y1_dB_lower = Y1_dB(:, 1:12*N+1);
% Y2_dB_lower = Y2_dB(:, 1:12*N+1);
% Y3_dB_lower = Y3_dB(:, 1:12*N+1);

% Normalisation
Y1_norm = Y1_dB_lower - max(max(Y1_dB_lower));
Y2_norm = Y2_dB_lower - max(max(Y2_dB_lower));
Y3_norm = Y3_dB_lower - max(max(Y3_dB_lower));

% compute Rmax
c = 299792458;
f_start = 2.405 * 10^9;
f_stop = 2.495 * 10^9;
delta_f = f_stop - f_start;
delta_R = c/(2*delta_f);
R_max = (N*delta_R)/2;
y =  linspace(0, R_max, half_size);

%plot the spectograms
imagesc(y, time, Y1_norm, [-50, 0])
title("Range with MS and no MTI")
xlabel("Range (m)")
ylabel("Time (s)")
xlim([0, 30])
colorbar
figure

time_2MTI = time(1:end-1);
imagesc(y, time_2MTI, Y2_norm, [-50, 0])
title("Range with MS and 2-pulse MTI")
xlabel("Range (m)")
ylabel("Time (s)")
xlim([0, 30])
colorbar
figure


time_3MTI = time(1:end-2);
imagesc(y, time_3MTI, Y3_norm, [-50, 0])
title("Range with MS and 3-pulse MTI")
xlabel("Range (m)")
ylabel("Time (s)")
xlim([0, 30])
colorbar

%%
min_dist = delta_R;      % Minimum spacing between peaks (m)
min_height = -25;   % Power threshold in dB

chosen_Y = Y3_norm;
chosen_time = time_3MTI;

r_est1 = zeros(size(chosen_Y,1),1); % Preallocate
r_est2 = zeros(size(chosen_Y,1),1);

for i = 1:size(chosen_Y,1)
    row = chosen_Y(i,:);
    [peakVals, locs] = findpeaks(row, y, ...
        'MinPeakHeight', min_height, ...
        'MinPeakDistance', min_dist);
    if isempty(peakVals)
        if i > 1
            r_est1(i) = r_est1(i-1);
            r_est2(i) = r_est2(i-1);
        else
            continue;
        end

    elseif isscalar(peakVals)
        r_est1(i) = locs(1);
        r_est2(i) = NaN;
    else
        % Take two strongest peaks
        [~, sortIdx] = maxk(peakVals, 2);
        peakLocs = locs(sortIdx);
        
        % Assign consistently: smaller range = Target 1, larger = Target 2
        r_sorted = sort(peakLocs);
        r_est1(i) = r_sorted(1);
        r_est2(i) = r_sorted(2);
    end
end
%%
figure;
plot(chosen_time, r_est1);
xlabel('Time (s)'); ylabel('Range (m)');
title('Range vs Time (Single Strongest Scatterer)');
xlim([0, 25]);
grid on;


figure;
plot(chosen_time, r_est2);
xlabel('Time (s)'); ylabel('Range (m)');
title('Range vs Time (Second Strongest Scatterer)');
xlim([0, 25]);
grid on;

figure;
plot(chosen_time, r_est1, 'DisplayName','Target 1 (close)'); hold on;
plot(chosen_time, r_est2, 'DisplayName','Target 2 (farther)');
xlabel('Time (s)');
ylabel('Range (m)');
title('Estimated Range of Two Targets');
legend show; grid on;

%%
half_size = round(size(Y3, 2)/2);
Y3_lower = Y3(:, 1:half_size);

% Plots without noise reduction
[val, idx] = max(Y3_lower, [], 2);
strongest_range = y(idx);
figure
plot(time_3MTI, val)
title("Size of scatter")
figure
plot(time_3MTI, strongest_range)
title("FMCW, 3MTI, 4N padding, MS clutter rejection")
legend('Strongest scatter')
xlabel('Time (s)')
ylabel('Range (m)')
%%

detection_threshold = 0.01;
max_vel = 5;

range = zeros(size(Y3_lower,1),1);
velocity = zeros(size(Y3_lower,1),1);

target_detected = false;
latest_update = 0;
for i = 1:size(Y3_lower, 1)
    row = Y3_lower(i, :);
    if target_detected
        reasonable_vel = false;
        [val, idx] = max(row);
        delta_length = y(idx) - range(latest_update);
        delta_time = time_3MTI(i) - time_3MTI(latest_update);
        vel = delta_length / delta_time;
        %disp(vel)
        if abs(vel) < max_vel % | i >= size(Y3_lower, 1)
            for j = latest_update:i
                velocity(j) = vel;
                range(j) = range(latest_update) + ...
                    vel * (time_3MTI(j) - time_3MTI(latest_update));
            end
            latest_update = i;
        end

    else
        if max(row) > detection_threshold
            target_detected = true;
            [val, idx] = max(row);
            range(i) = y(idx);
            latest_update = i;
            first_target = i;
        end
    end
end
range(1:first_target-1) = [];
velocity(1:first_target-1) = [];
hold on
plot(time_3MTI, range)
legend('Strongest scatter', 'Noise reduced target')
xlabel('Time (s)')
ylabel('Range (m)')
%%
figure
plot(time_3MTI, velocity)
hold on
velocity_smooth = movmean(velocity, 500);
plot(time_3MTI, velocity_smooth)
title("FMCW, 3MTI, 4N padding, MS clutter rejection")
xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend('Noise reduced velocity', 'Noise reduced velocity moving average')
