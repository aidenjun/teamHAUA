%% FINAL DEMO – Noise Signals / Audio Filtering (NSAF)
% Team HAUA
% Clean → Noisy → Analysis → Filtering → Comparison

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

%% 4. FFT & Spectrogram (Noise Identification)
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
title('FFT – Noisy'); grid on; xlim([0 Fs/2])

figure;
subplot(2,1,1)
spectrogram(clean, hamming(256), 128, 512, Fs, 'yaxis');
title('Clean Spectrogram'); ylim([0 8]);

subplot(2,1,2)
spectrogram(noisy, hamming(256), 128, 512, Fs, 'yaxis');
title('Noisy Spectrogram'); ylim([0 8]);

%% 5. Filter Design
notch = designfilt('bandstopiir','FilterOrder',4,...
    'HalfPowerFrequency1',58,'HalfPowerFrequency2',62,...
    'DesignMethod','butter','SampleRate',Fs);

lpf = designfilt('lowpassiir','FilterOrder',6,...
    'HalfPowerFrequency',4000,...
    'DesignMethod','butter','SampleRate',Fs);

%% 6. Time-Domain Filtering
clean_time = filtfilt(lpf, filtfilt(notch, noisy));

%% 7. Frequency-Domain Filtering
[b1,a1] = tf(notch);
[b2,a2] = tf(lpf);

delta = zeros(N_fft,1);
delta(1) = 1;

h = filter(b2, a2, filter(b1, a1, delta));
H = fft(h, N_fft);

X = fft(noisy, N_fft);
clean_freq = real(ifft(X .* H));
clean_freq = clean_freq(1:N);

%% 8. Spectrogram After Filtering
figure;
subplot(3,1,1)
spectrogram(noisy, hamming(256),128,512,Fs,'yaxis');
title('Noisy'); ylim([0 8]);

subplot(3,1,2)
spectrogram(clean_time, hamming(256),128,512,Fs,'yaxis');
title('Filtered (Time)'); ylim([0 8]);

subplot(3,1,3)
spectrogram(clean_freq, hamming(256),128,512,Fs,'yaxis');
title('Filtered (Freq)'); ylim([0 8]);

%% 9. Metrics
fprintf('SNR before: %.2f dB\n', snr(clean, noisy-clean));
fprintf('SNR after:  %.2f dB\n', snr(clean, clean_time-clean));
fprintf('Correlation: %.6f\n', corr(clean_time, clean_freq));

%% 10. Save Output Audio
audiowrite("clean_speech_filtered_time.wav", clean_time, Fs);
audiowrite("clean_speech_filtered_freq.wav", clean_freq, Fs);

disp('FINAL DEMO COMPLETE');
