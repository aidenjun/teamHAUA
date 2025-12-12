clear; clc; close all;

[clean, Fs] = audioread("wavfile.wav");  

if size(clean,2) > 1
    clean = mean(clean, 2);
end

clean = clean / max(abs(clean));

f_hum = 60;         
hum_scale = 0.1;    
[with_hum, hum] = add_hum(clean, Fs, f_hum, hum_scale);

snr_dB = 15;         
[with_hiss, hiss] = add_hiss(clean, snr_dB);

with_hum_and_hiss = with_hum + hiss;
Tplot = 0.05;                       
Nplot = min(length(clean), round(Tplot*Fs));
t = (0:Nplot-1)/Fs;

figure;
subplot(3,1,1);
plot(t, clean(1:Nplot));
xlabel('Time (s)');
ylabel('Amplitude');
title('Clean Audio (Time Domain)');
grid on;

subplot(3,1,2);
plot(t, with_hum_and_hiss(1:Nplot));
xlabel('Time (s)');
ylabel('Amplitude');
title(sprintf('Noisy Audio (Hum %d Hz + Hiss, SNR = %d dB)', f_hum, snr_dB));
grid on;

subplot(3,1,3);
plot(t, with_hum_and_hiss(1:Nplot) - clean(1:Nplot));
xlabel('Time (s)');
ylabel('Amplitude');
title('Noise Component (Hum + Hiss) in Time Domain');
grid on;


audiowrite("clean_speech_clean.wav", clean, Fs);
audiowrite("clean_speech_hum.wav", with_hum, Fs);
audiowrite("clean_speech_hiss.wav", with_hiss, Fs);
audiowrite("clean_speech_hum_hiss.wav", with_hum_and_hiss, Fs);

disp('audio files generated');
