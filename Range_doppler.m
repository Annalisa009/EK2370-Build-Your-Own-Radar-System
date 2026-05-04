clear all; close all; clc
% 0. Draw an empty plot to pre-load the plot tool
imagesc([], [], [], [-50, 0]);
xlabel("Range (m)");
ylabel("Velocity (m/s)");
xlim([0, 100]);
colorbar;
drawnow;

% 1. Read audio file
[data, Fs] = audioread("fmcw_range_doppler.wav");

% 2. Separate sync data and back scatter data
sync_data = data(:, 2);
back_scatter_data = data(:, 1);
%plot(sync_data)
%figure
%plot(back_scatter_data)

% 3. Place everything coming up in a while loop
idx = 1;
max_idx = size(sync_data, 1);
chunk = round(1 * Fs);

strongest_range = [];
strongest_vel = [];
counter = 0;
while idx + chunk <= max_idx
    counter = counter + 1;
    chunked_sync_data = sync_data(idx: idx + chunk);
    chunked_back_scatter_data = back_scatter_data(idx: idx + chunk);
    
    % 4. Binarize sync data
    threshold = 0;
    chunked_sync_data = chunked_sync_data > threshold;

    % 4 (again). Split the back scatter data
    Tp = 0.005;
    chirp_duration = round(Tp*Fs);
    combined_data = chunked_sync_data.*chunked_back_scatter_data;
    d = diff([0 transpose(chunked_sync_data) 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    num_segments = length(starts);
    seg_lengths = ends - starts + 1;
    max_length = max(seg_lengths);
    segments = zeros(num_segments, max_length);
    
    for k = 1:num_segments
        seg = combined_data(starts(k): min(ends(k), starts(k) + chirp_duration));
        segments(k, 1:length(seg)) = seg;
    end

    % 5. MS clutter rejection
    col_mean = mean(segments, 1);
    segments = segments - col_mean;
    
    % 6. 2 pulse MTI
    data_2MTI = segments(2:size(segments,1),:) - segments(1:size(segments,1)-1,:);
    %data_2MTI = segments;
    
    % 7. IFFT, no zero padding
    %N = chirp_duration;
    %N = half_size*2;
    Y1 = ifft(data_2MTI, [], 2);
    Y2 = ifft(Y1, [], 1);

    Y2_dB = 20 * log10(abs(Y2));
    
    % 8. Keep lower part of matrix
    half_size = floor(size(Y2_dB, 2)/2);
    N = half_size*2;
    Y2_dB_lower = Y2_dB(:, 1:half_size);

    % 9. Normalization
    Y2_norm = Y2_dB_lower - max(max(Y2_dB_lower));

    % 11 (no 10). Create data for y-axis
    c = 299792458;
    f_start = 2.405 * 10^9;
    f_stop = 2.495 * 10^9;
    f_c = (f_start + f_stop) / 2;
    f_D = linspace(0, 1/(2*Tp), size(Y2_norm, 1));
    velocity_array = (c/(2*f_c)) * f_D;

    % 12. Create data for x-axis
    delta_f = f_stop - f_start;
    delta_R = c/(2*delta_f);
    R_max = (N*delta_R)/2;
    range_array = linspace(0, R_max, size(Y2_norm, 2));

    % 13. Generate final spectrogram
    imagesc(range_array, velocity_array, Y2_norm, [-50, 0]);
    xlabel("Range (m)");
    ylabel("Velocity (m/s)");
    xlim([0, 40]);
    colorbar;
    drawnow;

    %
    Y2_lower = abs(Y2(:, 1:half_size));
    [val, linearIndex] = max(Y2_lower(:));
    [row, col] = ind2sub(size(Y2_lower), linearIndex);
    strongest_range = [strongest_range, range_array(col)];
    strongest_vel = [strongest_vel, velocity_array(row)];


    % 14. Increment the index
    idx = idx + chunk;


end
time_axis = linspace(0, counter*chunk/Fs, counter);
figure
plot(time_axis, strongest_vel)
title("FMCW, 2MTI, no zero padding, MS clutter rejection")
legend('Strongest scatter')
xlabel('Time(s)')
ylabel('Velocity(m/s)')
figure
plot(time_axis, strongest_range)
title("FMCW, 2MTI, no zero padding, MS clutter rejection")
legend('Strongest scatter')
xlabel('Time(s)')
ylabel('Range(m)')