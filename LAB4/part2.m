clc;
clear;

%% Load Data
load("E:\6th Semester\MISP Lab\MyLab\LAB4\SSVEP_EEG.mat");

% Extract data
eeg_data = SSVEP_Signal';  
trigger_indices = Event_samples;
fs = 256; 
num_channels = size(eeg_data, 2);  
channel_names = {'Pz', 'Oz', 'P7', 'P8', 'O2', 'O1'};  

fprintf('Data dimensions: %d samples × %d channels\n', size(eeg_data, 1), size(eeg_data, 2));
fprintf('Number of triggers: %d\n', length(trigger_indices));


%% Part (a): Bandpass Filtering (1-40 Hz) for ALL channels

% Design Butterworth bandpass filter
filter_order = 4;
low_cutoff = 1;   
high_cutoff = 40; 
[b, a] = butter(filter_order, [low_cutoff high_cutoff]/(fs/2), 'bandpass');

% Apply filter to ALL channels
filtered_data = zeros(size(SSVEP_Signal));
for ch = 1:num_channels
    filtered_data(ch, :) = filtfilt(b, a, SSVEP_Signal(ch, :));
end
figure('Name', 'Part (a): Original vs Filtered Signal - All Channels', 'Position', [50 50 1400 900]);
t = (0:size(SSVEP_Signal, 2)-1) / fs;

for ch = 1:num_channels
    subplot(3, 2, ch);
    hold on;
    plot(t, SSVEP_Signal(ch, :), 'b', 'LineWidth', 0.5, 'DisplayName', 'Original');
    plot(t, filtered_data(ch, :), 'r', 'LineWidth', 1.2, 'DisplayName', 'Filtered (1-40 Hz)');
    hold off;
    
    title(['Channel: ' channel_names{ch}]);
    xlabel('Time (seconds)');
    ylabel('Amplitude (μV)');
    legend('Location', 'best');
    grid on;
end
%% Part B) Extract trials (15 trials × 15 stimuli)
trial_length = 4 * fs;  % 4 seconds
num_trials = length(trigger_indices);
trials = zeros(trial_length, num_channels, num_trials);

for trial = 1:num_trials
    start_idx = trigger_indices(trial);
    end_idx = start_idx + trial_length - 1;
    
    if end_idx <= size(filtered_data, 2)  
        trials(:, :, trial) = filtered_data(:, start_idx:end_idx)';
    else
        warning('Trial %d exceeds data length. Skipping.', trial);
    end
end

fprintf('Extracted %d trials with length %d samples.\n', num_trials, trial_length);

% Plot example trials for all channels
figure('Name', 'Part (b): Example Trials - All Channels', 'Position', [50 50 1400 900]);
t_trial = (0:trial_length-1) / fs;

for ch = 1:num_channels
    subplot(3, 2, ch);
    hold on;
    for trial = 1:min(3, num_trials)  
        plot(t_trial, trials(:, ch, trial), 'DisplayName', sprintf('Trial %d', trial));
    end
    hold off;
    
    title(['Channel: ' channel_names{ch}]);
    xlabel('Time (seconds)');
    ylabel('Amplitude (μV)');
    legend('Location', 'best');
    grid on;
end


%% Part C) Calculate and plot PSD for each trial
pwelch_window = 2 * fs;  
overlap = pwelch_window / 2; 
nfft = 2^nextpow2(pwelch_window);

% Calculate PSD for each trial and each channel
psd_all = cell(num_trials, num_channels);
freq_all = [];

for trial = 1:num_trials
    for ch = 1:num_channels
        [psd, freq] = pwelch(trials(:, ch, trial), pwelch_window, overlap, nfft, fs);
        psd_all{trial, ch} = psd;
        if isempty(freq_all)
            freq_all = freq;
        end
    end
end

fprintf('PSD calculation completed.\n');

% Plot PSD for all trials - each channel in a separate figure
for ch = 1:num_channels
    figure('Name', sprintf('PSD All Trials - %s', channel_names{ch}));
    hold on;
    colors = lines(num_trials);
    
    for trial = 1:num_trials
        plot(freq_all, 10*log10(psd_all{trial, ch}), 'Color', colors(trial,:), ...
             'DisplayName', sprintf('Trial %d', trial));
    end
    
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title(sprintf('PSD All Trials - Channel %s', channel_names{ch}));
    xlim([0 50]);
    grid on;
    legend('Location', 'best', 'NumColumns', 3);
    hold off;
end

%% Average PSD across all trials
avg_psd = zeros(length(freq_all), num_channels);
for ch = 1:num_channels
    psd_matrix = cell2mat(cellfun(@(x) x, psd_all(:, ch), 'UniformOutput', false)');
    avg_psd(:, ch) = mean(psd_matrix, 2);
end

% Plot average PSD for all channels
figure('Name', 'Average PSD - All Channels');
num_rows = ceil(num_channels / 2);
for ch = 1:num_channels
    subplot(num_rows, 2, ch);
    plot(freq_all, 10*log10(avg_psd(:, ch)), 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title(sprintf('Average PSD - %s', channel_names{ch}));
    xlim([0 50]);
    grid on;
end

fprintf('Analysis completed.\n');
