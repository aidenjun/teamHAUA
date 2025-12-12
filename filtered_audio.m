clear; clc; close all;

[noisy, Fs] = audioread("clean_speech_hum_hiss.wav");
[clean, ~] = audioread("clean_speech_clean.wav");
noisy = noisy / max(abs(noisy));
clean = clean / max(abs(clean));

%Notch filter design (removes 60Hz hum)
notch_filter = designfilt('bandstopiir', 'FilterOrder',4,'HalfPowerFrequency1', 58,...
    'HalfPowerFrequency2', 62, 'DesignMethod', 'butter', 'SampleRate', Fs);


%Low pass filter design (removes high frequency hiss)
lpf = designfilt('lowpassiir','FilterOrder', 6,'HalfPowerFrequency', 4000, ...
    'DesignMethod', 'butter','SampleRate', Fs);


%Time Domain Convolution Application
    %Applyings filter using filtfilt()
tic;
noisy_notch = filtfilt(notch_filter, noisy);
noisy_clean_time = filtfilt(lpf, noisy_notch);
time_domain_duration = toc;
fprintf('Time Domain filtering duration (filtfilt): %.4f seconds\n', time_domain_duration);


%Frequency Domain Application
tic;
    % Get filter frequency responses at appropriate length
N = length(noisy);
[b_notch, a_notch] = tf(notch_filter);
[b_low, a_low] = tf(lpf);

impulse_response_length = max(length(b_notch), length(b_low)) * 100;

    % Use FFT length for linear convolution
N_fft = 2^nextpow2(N + impulse_response_length);

    % Create a delta function
delta = zeros(N_fft, 1);
delta(1) = 1;

    % Apply both filters to delta (forward filtering only, not filtfilt)
imp_response_notch = filter(b_notch, a_notch, delta);
imp_response_combined = filter(b_low, a_low, imp_response_notch);

    % Take FFT of the combined impulse response
H_combined_fft = fft(imp_response_combined, N_fft);

    % Compute FFT of noisy signal (with zero-padding for linear convolution)
X_noisy = fft(noisy, N_fft);

    % Apply filter in frequency domain (multiplication)
X_filtered = X_noisy .* H_combined_fft;

    % Convert back to time domain and keep only the valid part
noisy_clean_freq_full = real(ifft(X_filtered));
noisy_clean_freq = noisy_clean_freq_full(1:N);

freq_time = toc;
fprintf('Frequency Domain filtering time: %.4f seconds\n', freq_time);

%Time-domain comparison of all signals
figure('Position', [100, 100, 1200, 800]);

    % Plot all signals in time domain
t = (0:length(clean)-1)/Fs;

subplot(4,1,1);
plot(t, clean, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Clean Audio (Time Domain)');
grid on;

subplot(4,1,2);
plot(t, noisy, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Noisy Audio with 60Hz Hum + High-Frequency Hiss');
grid on;

subplot(4,1,3);
plot(t, noisy_clean_time, 'g', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Audio (Time-domain Method)');
grid on;

subplot(4,1,4);
plot(t, noisy_clean_freq, 'm', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Audio (Frequency-domain Method)');
grid on;

%FFT magnitude spectrum comparison
N_fft = 2^nextpow2(length(clean));
f = (0:N_fft-1)*(Fs/N_fft);
f = f(1:floor(N_fft/2));

    % Compute FFT magnitudes
fft_clean = abs(fft(clean, N_fft));
fft_noisy = abs(fft(noisy, N_fft));
fft_filtered_time = abs(fft(noisy_clean_time, N_fft));
fft_filtered_freq = abs(fft(noisy_clean_freq, N_fft));

    %Extract positive frequencies
mag_clean = fft_clean(1:floor(N_fft/2));
mag_noisy = fft_noisy(1:floor(N_fft/2));
mag_filtered_time = fft_filtered_time(1:floor(N_fft/2));
mag_filtered_freq = fft_filtered_freq(1:floor(N_fft/2));

    %Plot FFT comparison
figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(f, mag_clean, 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Magnitude: Original Clean Audio');
grid on;
xlim([0, Fs/2]);

subplot(2,2,2);
plot(f, mag_noisy, 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Magnitude: Noisy Audio');
grid on;
xlim([0, Fs/2]);

subplot(2,2,3);
plot(f, mag_filtered_time, 'g', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Magnitude: Filtered (Time-domain)');
grid on;
xlim([0, Fs/2]);

subplot(2,2,4);
plot(f, mag_filtered_freq, 'm', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Magnitude: Filtered (Frequency-domain)');
grid on;
xlim([0, Fs/2]);

%Combined comparison plots
figure('Position', [100, 100, 1200, 800]);

    %FFT magnitude comparison on single plot
subplot(2,1,1);
plot(f, mag_clean, 'b', 'LineWidth', 1.5, 'DisplayName', 'Original Clean');
hold on;
plot(f, mag_noisy, 'r', 'LineWidth', 1.5, 'DisplayName', 'Noisy');
plot(f, mag_filtered_time, 'g', 'LineWidth', 1.5, 'DisplayName', 'Filtered (Time)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Magnitude Comparison: All Signals');
legend('Location', 'best');
grid on;
xlim([0, Fs/2]);

% Difference between filtering methods
subplot(2,1,2);
plot(t, noisy_clean_time - noisy_clean_freq, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude Difference');
title('Difference Between Time-domain and Frequency-domain Filtering Methods');
grid on;

%Correlation Between filtering methods
correlation = corr(noisy_clean_time, noisy_clean_freq);
max_diff = max(abs(noisy_clean_time - noisy_clean_freq));
rms_diff = rms(noisy_clean_time - noisy_clean_freq);

fprintf('METHOD COMPARISON:\n');
fprintf('  • Correlation between methods: %.6f\n', correlation);
fprintf('  • Maximum absolute difference: %.6f\n', max_diff);
fprintf('  • RMS difference: %.6f\n', rms_diff);
fprintf('\n');

%Save filtered audio files
audiowrite('clean_speech_filtered_time.wav', noisy_clean_time, Fs);
audiowrite('clean_speech_filtered_freq.wav', noisy_clean_freq, Fs);