clc; clear; close all;

%% Parameters
Fs = 44100;              % Sampling frequency (match radar output!)
nBits = 16;              % Bits per sample
nChannels = 2;           % Stereo input (radar on channel 1)
Tp = 0.5;                % Sweep time [s]
N = round(Tp*Fs);        % Samples per sweep
c = 3e8;                 % Speed of light
fc = 2.445e9;            % Radar carrier frequency

%% 
audiodevinfo
%% Create recorder object
recObj = audiorecorder(Fs, 16, 1, DeviceID);
disp('>> Starting real-time radar acquisition...');
record(recObj);          % start recording in background

%% Velocity axis
f_axis = linspace(0, Fs/2, 2*N);
v_axis = (c/(2*fc)) * f_axis;

%% Initialize spectrogram buffer
maxSweeps = 200;   % how many sweeps to display in the live spectrogram
Y_buffer = nan(maxSweeps, 2*N);  % store FFT magnitudes
t_axis = (0:maxSweeps-1) * Tp;   % time axis

%% Prepare plot
figure;
hImg = imagesc(v_axis, t_axis, Y_buffer, [-55 0]); % dB range
set(gca, 'YDir', 'normal');
xlabel('Velocity (m/s)');
ylabel('Time (s)');
title('Real-time Spectrogram');
colormap(jet);
colorbar;
ylim([0 maxSweeps*Tp]);
xlim([0 30]);

%% Real-time loop
sweepCount = 0;
while isrecording(recObj)
    pause(Tp);   % wait one sweep duration

    % Fetch all recorded samples so far
    data = getaudiodata(recObj);
    radar = data(:,1);   % use first channel

    if length(radar) >= N
        % Take only the last sweep
        sweep = radar(end-N+1:end);

        % Mean subtraction
        sweep_ms = sweep - mean(sweep);

        % FFT
        Y = fft(sweep_ms, 4*N);
        Y_dB = 20*log10(abs(Y(1:2*N)));
        Y_dB = Y_dB - max(Y_dB); % normalize per sweep

        % Update spectrogram buffer
        sweepCount = sweepCount + 1;
        rowIdx = mod(sweepCount-1, maxSweeps) + 1; % circular buffer
        Y_buffer(rowIdx,:) = Y_dB;

        % Update plot
        set(hImg, 'CData', Y_buffer, ...
                  'YData', (0:maxSweeps-1)*Tp, ...
                  'XData', v_axis);
        ylim([0 maxSweeps*Tp]); % sliding time window
        drawnow;
    end
end
