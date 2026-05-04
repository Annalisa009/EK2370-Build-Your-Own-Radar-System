Fs = 44100;              % Adjust to your soundcard sampling rate
Tp = 0.005;              % Chirp duration (s)
chirp_samples = round(Tp * Fs);

%Run this command to see which audio devices are connected to see the USB sound card
%audiodevinfo

% Setup audio recorder (stereo: ch1=backscatter, ch2=sync)
deviceID = 1;  % (example, replace with your USB sound card’s ID from step 1)
recObj = audiorecorder(Fs, 16, 2, deviceID);
record(recObj);

% Radar parameters
c = 299792458;
f_start = 2.405e9;
f_stop  = 2.495e9;
delta_f = f_stop - f_start;
delta_R = c/(2*delta_f);
R_max   = (chirp_samples * delta_R)/2;
y = linspace(0, R_max, 2*chirp_samples);

% Spectrogram buffer
max_frames = 200; % number of chirps shown in plot
spec_buffer = -50*ones(max_frames, 2*chirp_samples); % initialize with -50 dB

figure;
h = imagesc(y, (0:max_frames-1)*Tp, spec_buffer, [-50, 0]);
xlabel("Range (m)");
ylabel("Time (s)");
title("Real-time FMCW Range Spectrogram (2-pulse MTI)");
colorbar;
xlim([0,100]);

frame_idx = 0;
prev_chirp = []; % store previous chirp for MTI

%% Continuous processing loop
while true
    pause(Tp);  % wait approx. one chirp
    
    % Get available audio
    data = getaudiodata(recObj);
    if size(data,1) < chirp_samples
        continue; % not enough samples yet
    end
    
    % Separate channels
    backscatter = data(:,1);
    sync_data   = data(:,2) > 0; % threshold
    
    % Detect chirp start/stop using sync
    d = diff([0; sync_data; 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    
    if isempty(starts)
        continue; % no chirp detected yet
    end
    
    % Process the last complete chirp
    k = starts(end);
    if k+chirp_samples-1 > length(backscatter)
        continue; % wait for full chirp
    end
    chirp_segment = backscatter(k:k+chirp_samples-1);
    
    % Mean subtraction (MS clutter rejection)
    chirp_segment = chirp_segment - mean(chirp_segment);
    
    % Apply 2-pulse MTI if previous chirp exists
    if ~isempty(prev_chirp)
        mti_chirp = chirp_segment - prev_chirp;
        
        % FFT and normalization
        Y = abs(ifft(mti_chirp, 4*chirp_samples));
        Y_dB = 20*log10(Y(1:2*chirp_samples+1));
        Y_norm = Y_dB - max(Y_dB);
        
        % Update buffer
        frame_idx = frame_idx + 1;
        if frame_idx > max_frames
            spec_buffer = circshift(spec_buffer, -1, 1);
            spec_buffer(end,:) = Y_norm;
        else
            spec_buffer(frame_idx,:) = Y_norm;
        end
        
        % Update spectrogram plot
        set(h, 'CData', spec_buffer);
        drawnow;
    end
    
    % Store current chirp as previous
    prev_chirp = chirp_segment;
end
