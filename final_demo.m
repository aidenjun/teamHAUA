%% FINAL DEMO – Noise Signals / Audio Filtering (NSAF)
% Team HAUA
% Base MATLAB ONLY (no toolboxes)

clear; clc; close all;

%% 1. Load Clean Audio
[clean, Fs] = audioread("wavfile.wav");
if size(clean,2) > 1
    clean = mean(clean,2);
end
clean = clean / max(abs(clean));
t = (0:length(clean)-1)/Fs;

fprintf('Loaded audio: %.2f s, Fs = %d Hz\n', length(clean)/Fs, Fs);

%% 2. Add Noise (60 Hz Hum + Hiss)
f_hum = 60;
hum_scale = 0.1;
snr_dB = 15;

[with_hum, ~] = add_sin_hum(clean, Fs, f_hum, hum_scale);
[~, hiss] = add_highfreq_hiss(clean, snr_dB);

noisy = with_hum + hiss;
noisy = noisy / max(abs(noisy));

audiowrite("clean_speech_clean.wav", clean, Fs);
audiowrite("clean_speech_hum_hiss.wav", noisy, Fs);

%% 3. Time-Domain Comparison (Zoomed)
Tplot = 0.05;
Nplot = min(length(clean), round(Tplot*Fs));

figure;
subplot(3,1,1)
plot(t(1:Nplot), clean(1:Nplot))
title('Clean Audio'); grid on;

subplot(3,1,2)
plot(t(1:Nplot), noisy(1:Nplot))
title('Noisy Audio (Hum + Hiss)'); grid on;

subplot(3,1,3)
plot(t(1:Nplot), noisy(1:Nplot) - clean(1:Nplot))
title('Noise Component'); grid on;

%% 4. FFT Analysis (Noise Identification)
N = length(clean);
N_fft = 2^nextpow2(N);

f = (0:N_fft-1)*(Fs/N_fft);
f = f(1:N_fft/2);

FFT_clean = abs(fft(clean, N_fft));
FFT_noisy = abs(fft(noisy, N_fft));

figure;
subplot(2,1,1)
plot(f, FFT_clean(1:N_fft/2))
title('FFT – Clean'); grid on; xlim([0 Fs/2])

subplot(2,1,2)
plot(f, FFT_noisy(1:N_fft/2))
title('FFT – Noisy (60 Hz Spike + Hiss)')
grid on; xlim([0 Fs/2])

%% 5. Frequency-Domain Filtering (Toolbox-Free)
X = fft(noisy, N_fft);
H = ones(N_fft,1);

% --- 60 Hz notch (manual) ---
bw = 2; % Hz bandwidth
idx_60 = round(f_hum * N_fft / Fs);
bw_bins = round(bw * N_fft / Fs);

H(idx_60-bw_bins : idx_60+bw_bins) = 0;
H(N_fft-(idx_60+bw_bins) : N_fft-(idx_60-bw_bins)) = 0;

% --- Low-pass filter at 4 kHz ---
f_cut = 4000;
idx_cut = round(f_cut * N_fft / Fs);

H(idx_cut : end-idx_cut) = 0;

% Apply filter
Y = X .* H;
clean_freq = real(ifft(Y));
clean_freq = clean_freq(1:N);

%% 6. FFT Comparison After Filtering
FFT_filtered = abs(fft(clean_freq, N_fft));

figure;
plot(f, FFT_noisy(1:N_fft/2),'r'); hold on;
plot(f, FFT_filtered(1:N_fft/2),'g');
legend('Noisy','Filtered');
title('FFT Comparison After Filtering');
grid on; xlim([0 Fs/2])

%% 7. Manual SNR Calculation (Base MATLAB)
signal_power = mean(clean.^2);

noise_before = noisy - clean;
noise_after  = clean_freq - clean;

noise_power_before = mean(noise_before.^2);
noise_power_after  = mean(noise_after.^2);

snr_before = 10*log10(signal_power / noise_power_before);
snr_after  = 10*log10(signal_power / noise_power_after);

fprintf('SNR before filtering: %.2f dB\n', snr_before);
fprintf('SNR after filtering:  %.2f dB\n', snr_after);

%% 8. Save Output Audio
audiowrite("clean_speech_filtered_freq.wav", clean_freq, Fs);

disp('FINAL DEMO COMPLETE');

