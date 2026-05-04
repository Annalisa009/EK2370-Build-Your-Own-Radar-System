clear; clc; close all

%% Load data and set parameters
use_IQ = true;

use_LowIF = true;

if use_IQ
    [I_data, Fs] = audioread('SAR_data/breath_lowIF_I.wav');
    [Q_data, Fs] = audioread('SAR_data/breath_lowIF_Q.wav');
else
    [data, Fs] = audioread('SAR_data/single_lowIF_I.wav');
end

c  = 299792458; % Speed of light
f0 = 5.8e9; % Base radar frequency

if use_LowIF
    f0 = f0 + 4e3;
end

if use_IQ
    data = complex(I_data,Q_data);
end

% Remove first 3 seconds
data = data(3e6:end);
if use_IQ
    I_data = I_data(3e6:end);
    Q_data = Q_data(3e6:end);
end

%% Plot the recorded data
% t = (0:length(data)-1)/Fs;
% figure
% if use_IQ
%     subplot(2,1,1); plot(t, I_data); grid on;
%     xlabel('Time (s)'); ylabel('Amplitude'); title('Recorded I_data')
%     subplot(2,1,2); plot(t, Q_data); grid on;
%     xlabel('Time (s)'); ylabel('Amplitude'); title('Recorded Q_data')
% else
%     plot(t, data); grid on;
%     xlabel('Time (s)'); ylabel('Amplitude'); title('Recorded data')
% end

%% Data processing
if ~use_IQ
    data = hilbert(data);
end
phi = unwrap(angle(data));
R = (c/(4*pi*f0)) * phi;

%% Plot the 
t = (0:length(data)-1)/Fs;
figure
plot(t, 1e3*R);
grid on;
xlabel('Time (s)'); ylabel('Range (mm)');
title('Range')
exportgraphics(gcf,'myplot.pdf','ContentType','vector');

%%
N = length(data);
F = fft(R, 1*N);
F = F(1:0.5*N);
f = linspace(0, Fs/2, 0.5*N); 
figure
plot(f, 20*log10(abs(F)), 'LineWidth', 2)
xlim([0,3])
grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('FFT of range')