clc; clear all; close all;

%% Step 0: Read audio file
use_IQ = false;

for gain = [60, 70]
    for use_blocker = [true, false]
        name = "With_PA/PA_gain" + gain;
        if use_blocker
            I_name = name + "_blocker_I.wav";
            Q_name = name + "_blocker_Q.wav";
        else
            I_name = name + "_I.wav";
            Q_name = name + "_Q.wav";
        end
        
        if use_IQ
            [I_data, Fs] = audioread(I_name);
            [Q_data, Fs] = audioread(Q_name);
        else
            disp(I_name)
            [data, Fs] = audioread(I_name);
        end
        
        %% Step 1: Inspect both channels
        % figure; hold on;
        % plot(I_data(:,1),'b');    % channel 1 (inverted for visibility)
        % plot(Q_data(:,1),'r');     % channel 2
        % xlabel('Sample number');
        % ylabel('Amplitude');
        % legend('Channel 1','Channel 2');
        % title('Comparison of both channels');
        
        %% Pick radar channel
        if use_IQ
            radar = complex(I_data(:,1), Q_data(:,1));
        else
            radar = data(:,1);          % Channel 1 contains the radar data
        end
        t = (0:length(radar)-1)/Fs;
        
        %% Split radar data into sweeps
        Tp = 0.05;                    % in seconds - pulse length 
        N = round(Tp * Fs);          % samples per sweep
        M = floor(length(radar)/N)-1;  % number of complete sweeps
        
        radar = radar(1:M*N);  % truncate radar to integer number of sweeps
        sweeps = reshape(radar, [N, M]).';   % size M x N
        
        fprintf('Matrix size: %d sweeps x %d samples\n', M, N);
        
        % %% Optional: plot first sweep
        % t_sweep = (0:N-1)/Fs;
        % figure;
        % plot(t_sweep, sweeps(1,:));
        % xlabel('Time (s)');
        % ylabel('Amplitude');
        % title('First radar sweep');
        
        %% MS clutter rejection 
        disp(mean(sweeps(:)))
        mean_val = mean(sweeps(:));   % global mean
        sweeps_ms = sweeps - mean_val;
        
        w = hann(N).';
        sweeps_ms = sweeps_ms .* w;
        
        % % Optional, just for visualization - Plot first sweep after MS
        % figure;
        % plot(t_sweep, sweeps_ms(1,:));
        % xlabel('Time (s)');
        % ylabel('Amplitude');
        % title('First sweep after MS clutter rejection (global mean)');
        
        %% FFT with zero-padding and log scale
        Y = fft(sweeps_ms, 4*N, 2);        % FFT across columns, zero-padding to 4*N
        if use_IQ
            Y = fftshift(Y, 2);
        end
        Y_dB = 20*log10(abs(Y));           % conversion to dB
        
        if ~use_IQ
            Y_dB_half = Y_dB(:, 1:2*N); % Keep up to Fs/2
        else
            Y_dB_half = Y_dB;
        end
        
        %% Maximum representable frequency
        Fmax = Fs/2;  % Nyquist
        
        %% Normalizations
        
        Y_norm1 = Y_dB_half - max(max(Y_dB_half)); % subtract maximum of entire lower-half matrix
        
        Y_norm2 = Y_dB_half - max(Y_dB_half,[],2);  % row-wise max subtraction
        
        %% Doppler frequency array to velocity array
        if use_IQ
            f_axis = linspace(-Fmax, Fmax, 4*N);
        else
            f_axis = linspace(0, Fmax, 2*N);
        end
        
        % CW radar conversion: v = c/(2*fc) * fD
        c = 3e8;        
        fc = 5.8e9;    
        v_axis = (c/(2*fc)) * f_axis;  
        
        %% Time axis for sweeps
        t_axis = linspace(0, Tp*M, M);  
        
        %% Step 9: Plot spectrograms
        
        % figure;
        % imagesc(v_axis, t_axis, Y_norm1);
        % set(gca, 'YDir', 'reverse');
        % axis xy;
        % xlabel('Velocity (m/s)');
        % ylabel('Time (s)');
        % title('Spectrogram: normalized by global max');
        % colorbar;
        % caxis([-45 0]);  
        % colormap(jet);
        % if use_IQ
        %     xlim([-8 8]);
        % else
        %     xlim([0 8]);
        % end
        
        figure;
        imagesc(v_axis, t_axis, Y_norm2);
        set(gca, 'YDir', 'reverse');
        axis xy;
        xlabel('Velocity (m/s)');
        ylabel('Time (s)');
        title_name = "Gain: " + gain;
        if use_blocker
            title_name = title_name + ", with DC blocker";
        else
            title_name = title_name + ", without DC blocker";
        end
        title(title_name);
        colorbar;
        caxis([-10 0]);   % set color scale from -55 dB (blue) to 0 dB (yellow/red)
        colormap(jet);
        if use_IQ
            xlim([-8 8]);
        else
            xlim([0 8]);
        end
        drawnow;
    end
end