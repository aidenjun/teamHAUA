clear; clc; close all;

% Clean and noisy audio files
[clean, Fs] = audioread("clean_speech_clean.wav");
[noisy, ~] = audioread("clean_speech_hum_hiss.wav");

% Ensure mono
if size(clean,2) > 1
    clean = mean(clean, 2);
end
if size(noisy,2) > 1
    noisy = mean(noisy, 2);
end

% FFT Analysis
N = length(clean);
f = (0:N-1) * (Fs/N);  % Frequency vector
f = f(1:floor(N/2));   % Only positive frequencies

% Compute FFT
fft_clean = fft(clean);
fft_noisy = fft(noisy);

% Magnitude spectrum (only positive frequencies)
mag_clean = abs(fft_clean(1:floor(N/2)));
mag_noisy = abs(fft_noisy(1:floor(N/2)));

% Plot FFT Comparison
figure('Position', [100, 100, 1200, 800]);

% Full spectrum view
subplot(3,2,1);
plot(f, mag_clean, 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Clean Audio - FFT Magnitude');
grid on;
xlim([0 Fs/2]);

subplot(3,2,2);
plot(f, mag_noisy, 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Noisy Audio - FFT Magnitude');
grid on;
xlim([0 Fs/2]);

% Zoomed view: 0-200 Hz to see 60 Hz spike
subplot(3,2,3);
plot(f, mag_clean, 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Clean Audio - Low Frequency Range');
grid on;
xlim([0 200]);

subplot(3,2,4);
plot(f, mag_noisy, 'r', 'LineWidth', 1.5);
hold on;
% Mark the 60 Hz spike
[~, idx_60] = min(abs(f - 60));
plot(f(idx_60), mag_noisy(idx_60), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
text(f(idx_60)+5, mag_noisy(idx_60), '60 Hz Hum', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Noisy Audio - 60 Hz Spike Visible');
grid on;
xlim([0 200]);

% Comparison plot
subplot(3,2,5);
plot(f, mag_clean, 'b', 'LineWidth', 1.5, 'DisplayName', 'Clean');
hold on;
plot(f, mag_noisy, 'r', 'LineWidth', 1.5, 'DisplayName', 'Noisy');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Comparison: Clean vs Noisy');
legend('Location', 'northeast');
grid on;
xlim([0 5000]);

% Difference (Noise component)
subplot(3,2,6);
plot(f, mag_noisy - mag_clean, 'k', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude Difference');
title('Noise Component in Frequency Domain');
grid on;
xlim([0 5000]);

% Spectrogram Analysis
figure('Position', [100, 100, 1200, 600]);

% Clean spectrogram
subplot(2,1,1);
spectrogram(clean, hamming(256), 128, 512, Fs, 'yaxis');
title('Clean Audio - Spectrogram');
colorbar;
ylim([0 8]);  % Focus on 0-8 kHz

% Noisy spectrogram
subplot(2,1,2);
spectrogram(noisy, hamming(256), 128, 512, Fs, 'yaxis');
title('Noisy Audio - Spectrogram (60 Hz Hum + High-Frequency Hiss Visible)');
colorbar;
ylim([0 8]);  % Focus on 0-8 kHz

% Summary statistics ouput
disp('=== Frequency Analysis Summary ===');
fprintf('Sampling Rate: %d Hz\n', Fs);
fprintf('Signal Length: %.2f seconds\n', N/Fs);
fprintf('Frequency Resolution: %.2f Hz\n', Fs/N);
fprintf('\n60 Hz Component:\n');
fprintf('  Clean magnitude at 60 Hz: %.4f\n', mag_clean(idx_60));
fprintf('  Noisy magnitude at 60 Hz: %.4f\n', mag_noisy(idx_60));
fprintf('  Increase due to hum: %.2fx\n', mag_noisy(idx_60)/mag_clean(idx_60));

disp('Frequency analysis complete - noise identified.');
