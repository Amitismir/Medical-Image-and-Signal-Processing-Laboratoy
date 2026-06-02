clear; clc; close all;

%% 1. Load data
data=load("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\MISP Lab - Source\Lab 5\Lab 5_data\normal.mat");

fprintf('Variables in the .mat file:\n');
vars = whos;
for i = 1:length(vars)
    fprintf('  %s (size: %s, class: %s)\n', vars(i).name, mat2str(vars(i).size), vars(i).class);
end
if exist('data', 'var') && isstruct(data)
    fprintf('\n"data" is a structure. Fields:\n');
    fields = fieldnames(data);
    for i = 1:length(fields)
        fprintf('  .%s (size: %s)\n', fields{i}, mat2str(size(data.(fields{i}))));
    end
        if isfield(data, 'normal')
        data_matrix = data.normal;  % This is the [75100×2 single] matrix
        fprintf('\nUsing data.normal (size: %s)\n', mat2str(size(data_matrix)));
    else
        for i = 1:length(fields)
            if isnumeric(data.(fields{i}))
                data_matrix = data.(fields{i});
                fprintf('\nUsing data.%s\n', fields{i});
                break;
            end
        end
    end
else
    error('No structure variable named "data" found in the file');
end

data_matrix = double(data_matrix);

time = data_matrix(:, 1);     % time vector (seconds)
ecg_raw = data_matrix(:, 2);  % raw ECG signal (volts)

Fs = 250;             % sampling frequency (Hz)
Ts = 1/Fs;            % sampling period

fprintf('\nData loaded successfully!\n');
fprintf('Matrix size: %d rows x %d columns\n', size(data_matrix, 1), size(data_matrix, 2));
fprintf('Total samples: %d\n', length(time));
fprintf('Time duration: %.2f seconds\n', time(end));
fprintf('Sampling frequency: %d Hz\n', Fs);

% Select clean (5-10 s) and noisy (last 60 s, pick 5-10 s within it)
idx_clean = (time >= 5) & (time <= 10);

if sum(idx_clean) == 0
    fprintf('\nWarning: No data between 5-10 seconds.\n');
    fprintf('Time range of data: [%.2f, %.2f] seconds\n', time(1), time(end));
    idx_clean = (time >= time(1)) & (time <= time(1) + 5);
    fprintf('Using first 5 seconds as clean segment instead.\n');
end

time_clean = time(idx_clean);
ecg_clean = ecg_raw(idx_clean);
t_end = time(end);
t_start_noisy = t_end - 60;   % beginning of noisy period

if t_start_noisy < time(1)
    t_start_noisy = time(1);
end

% Pick 5-10 seconds from the noisy period
idx_noisy = (time >= t_start_noisy + 5) & (time <= t_start_noisy + 10);

if sum(idx_noisy) == 0
    idx_noisy = (time >= time(end)-10) & (time <= time(end));
    fprintf('Using last 10 seconds as noisy segment.\n');
end

time_noisy = time(idx_noisy);
ecg_noisy = ecg_raw(idx_noisy);

fprintf('\nSegment info:\n');
fprintf('  Clean segment: %.2f to %.2f seconds (%d samples)\n', time_clean(1), time_clean(end), length(ecg_clean));
fprintf('  Noisy segment: %.2f to %.2f seconds (%d samples)\n', time_noisy(1), time_noisy(end), length(ecg_noisy));

% Frequency analysis using pwelch (PSD in dB)
figure('Name', 'Power Spectral Density (PSD)', 'Position', [100, 100, 900, 600]);

% Clean segment
[psd_clean, f_clean] = pwelch(ecg_clean, hamming(256), 128, 1024, Fs);
psd_clean_dB = 10*log10(psd_clean + eps);  % eps to avoid log10(0)

subplot(2,1,1);
plot(f_clean, psd_clean_dB, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('PSD of Clean ECG Segment (5-10 s)');
xlim([0 125]);
ylim([-80 20]);

% Noisy segment
[psd_noisy, f_noisy] = pwelch(ecg_noisy, hamming(256), 128, 1024, Fs);
psd_noisy_dB = 10*log10(psd_noisy + eps);

subplot(2,1,2);
plot(f_noisy, psd_noisy_dB, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('PSD of Noisy ECG Segment (5-10 s within last minute)');
xlim([0 125]);
ylim([-80 20]);

sgtitle('Frequency Content Comparison: Clean vs Noisy ECG');

%% 2. Design bandpass filter

low_cutoff = 0.5;   
high_cutoff = 50;  
order = 4;          % 4th order Butterworth

% Normalized frequencies (Nyquist = Fs/2 = 125 Hz)
Wn = [low_cutoff high_cutoff] / (Fs/2);

% Design bandpass filter
[b, a] = butter(order, Wn, 'bandpass');

% Frequency response
[H, w] = freqz(b, a, 1024, Fs);
mag_dB = 20*log10(abs(H) + eps);
phase_deg = unwrap(angle(H)) * 180/pi;

% Impulse response
impulse = [1, zeros(1, 200)];
impulse_response = filter(b, a, impulse);
t_impulse = (0:length(impulse_response)-1) / Fs;

% Plot filter characteristics
figure('Name', 'Filter Design', 'Position', [100, 100, 1000, 800]);

subplot(2,2,1);
stem(t_impulse, impulse_response, 'filled', 'MarkerSize', 3);
grid on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Impulse Response of Bandpass Filter');

subplot(2,2,2);
plot(w, mag_dB, 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Response (Magnitude)');
xlim([0 80]);
ylim([-60 5]);
hold on;
plot([low_cutoff low_cutoff], [-60 5], 'r--', 'LineWidth', 1);
plot([high_cutoff high_cutoff], [-60 5], 'r--', 'LineWidth', 1);
legend('Filter Response', 'Cutoff Frequencies', 'Location', 'southwest');

subplot(2,2,3);
plot(w, phase_deg, 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Frequency Response (Phase)');
xlim([0 80]);

subplot(2,2,4);
zplane(b, a);
title('Pole-Zero Plot');

sgtitle(sprintf('Bandpass Filter Characteristics (%.1f - %.0f Hz, %dth order Butterworth)', low_cutoff, high_cutoff, order));

%% 6. Filter the clean and noisy segments (using filtfilt for zero-phase filtering)
ecg_clean_filtered = filtfilt(b, a, ecg_clean);
ecg_noisy_filtered = filtfilt(b, a, ecg_noisy);

%% 7. Demonstrate filter performance
figure('Name', 'Filtering Results', 'Position', [100, 100, 900, 600]);

% Clean segment: before and after
subplot(2,1,1);
plot(time_clean, ecg_clean, 'b', 'LineWidth', 1); hold on;
plot(time_clean, ecg_clean_filtered, 'r', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Clean ECG Segment: Original vs Filtered');
legend('Original', 'Filtered (0.5-40 Hz)');
xlim([time_clean(1) time_clean(end)]);

% Noisy segment: before and after
subplot(2,1,2);
plot(time_noisy, ecg_noisy, 'b', 'LineWidth', 1); hold on;
plot(time_noisy, ecg_noisy_filtered, 'r', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Noisy ECG Segment: Original vs Filtered');
legend('Original (with muscle noise)', 'Filtered');
xlim([time_noisy(1) time_noisy(end)]);

sgtitle('Effect of Bandpass Filtering on ECG Signal');

%% 8. Quantitative evaluation of noise reduction
rms_noisy_before = rms(ecg_noisy);
rms_noisy_after = rms(ecg_noisy_filtered);
noise_reduction_dB = 20*log10(rms_noisy_before / rms_noisy_after);

rms_clean_before = rms(ecg_clean);
rms_clean_after = rms(ecg_clean_filtered);
signal_distortion_dB = 20*log10(rms_clean_after / rms_clean_before);

fprintf('\n===== Filter Performance =====\n');
fprintf('Noisy segment:\n');
fprintf('  RMS before filtering: %.4f V\n', rms_noisy_before);
fprintf('  RMS after filtering : %.4f V\n', rms_noisy_after);
fprintf('  Noise reduction     : %.2f dB\n', noise_reduction_dB);
fprintf('\nClean segment (distortion check):\n');
fprintf('  RMS before filtering: %.4f V\n', rms_clean_before);
fprintf('  RMS after filtering : %.4f V\n', rms_clean_after);
fprintf('  Signal distortion   : %.2f dB (should be small)\n', signal_distortion_dB);

%% 9. Additional plot: PSD of noisy segment before and after filtering
figure('Name', 'PSD Before and After Filtering', 'Position', [100, 100, 900, 600]);
[psd_noisy_before, f_psd] = pwelch(ecg_noisy, hamming(256), 128, 1024, Fs);
[psd_noisy_after, ~] = pwelch(ecg_noisy_filtered, hamming(256), 128, 1024, Fs);

plot(f_psd, 10*log10(psd_noisy_before + eps), 'b', 'LineWidth', 1.5); hold on;
plot(f_psd, 10*log10(psd_noisy_after + eps), 'r', 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Noisy ECG: PSD Before and After Bandpass Filtering');
legend('Before filtering', 'After filtering');
xlim([0 125]);
ylim([-80 20]);

% Add vertical lines for cutoff frequencies
hold on;
plot([low_cutoff low_cutoff], [-80 20], 'k--', 'LineWidth', 1);
plot([high_cutoff high_cutoff], [-80 20], 'k--', 'LineWidth', 1);
legend('Before filtering', 'After filtering', 'Cutoff frequencies');

%% 10. Display filter coefficients
fprintf('\nFilter coefficients (b - numerator):\n');
disp(b);
fprintf('Filter coefficients (a - denominator):\n');
disp(a);