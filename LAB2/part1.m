clc;
clear;
close all;
%% Loading the Lab2_1 data
addpath("E:\6th Semester\MISP Lab\MyLab\LAB2\Lab 2_data\Lab 2_data\Lab2_1")
%% 1.1 Plotting the Channels in Time-Domain
load('Lab2_1/X_org.mat');
load('Lab2_1/Electrodes.mat');
electrode_name = Electrodes.labels;
offset = max(abs(X_org(:)));
fs = 256;
disp_eeg(X_org, offset, fs, electrode_name);
title('Time-Domain Represeantation of All EEG Signals', Interpreter = 'latex');
%% 1.2 Plotting noise
load('Lab2_1/X_noise.mat');
load('Lab2_1/Electrodes.mat');
electrode_name = Electrodes.labels;
offset = max(abs(X_org(:)));
fs = 256;
disp_eeg(X_noise, offset, fs, electrode_name);
title('Time-Domain Represeantation of All Noise Signals on each EEG Channel', Interpreter = 'latex');
%% 1.3 Adding the Noise and Comparison

snr = [-5, -15];
X_observed = cell(1, length(snr));
% We use the average power in our calculations
P_signal = mean(X_org(:).^2); 
P_noise_raw = mean(X_noise(:).^2);
% Adding the noise is as follows: X_obs = X_org + alpha * X_noise where
% alpha is sqrt(P_signal/ (P_noise * 10^(SNR/10))
for i=1:length(snr)
    snr_db = snr(i);
    P_noise_desired = P_signal / (10^(snr_db/10));
    alpha = sqrt(P_noise_desired/P_noise_raw);
    % We should do the scaling here our raw noise already has a power
    % itself
    X_observed{i} = X_org + alpha * X_noise;
end
X_observed_snr5 = X_observed{1};
X_observed_snr15 = X_observed{2};

% Plotting the final noisy eeg signals
offset5 = max(abs(X_observed_snr5(:)));
offset15 = max(abs(X_observed_snr15(:)));
fs = 256;
disp_eeg(X_observed_snr5, offset5,fs,electrode_name);
title('EEG Signal with -5dB Noise', Interpreter = 'latex');
disp_eeg(X_observed_snr15,offset15,fs,electrode_name);
title('EEG Signal with -15dB Noise', Interpreter = 'latex');
%% 1.4
[F5,W5,K5] = COM2R(X_observed_snr5, 32);
[F15,W15,K15] = COM2R(X_observed_snr15, 32);

independent_source_signal5 = W5 * X_observed_snr5;
independent_source_signal15 = W15 * X_observed_snr15;
offset_5 = max(abs(independent_source_signal5(:)));
offset_15 = max(abs(independent_source_signal15(:)));
fs = 256;
disp_eeg(independent_source_signal5, offset_5, fs, electrode_name);
title('Independent Components with -5 dB SNR', Interpreter = 'latex');
disp_eeg(independent_source_signal15,offset_15, fs, electrode_name);
title('Independent Components with -15 dB SNR', Interpreter = 'latex');
%% 1.5 and 1.6

spikey_components_5 = [2, 5, 6, 9, 11, 12, 18];
selected_5 = independent_source_signal5(spikey_components_5, :);
reconstructed_5 = F5(:,spikey_components_5) * selected_5;
spikey_components_15 = [ 10, 20, 30];
selected_15 = independent_source_signal15(spikey_components_15, :);
reconstructed_15 = F15(:,spikey_components_15) * selected_15;

offset5 = max(abs(reconstructed_5(:)));
offset15 = max(abs(reconstructed_15(:)));
disp_eeg(reconstructed_5, offset5, fs, electrode_name);
title('The Reconstructed signal from -5 dB Noise');
disp_eeg(reconstructed_15, offset15, fs, electrode_name);
title('The Reconstructed signal from -15 dB Noise');
%% 1.7

ch13_5 = [reconstructed_5(13, :), X_org(13, :), X_observed_snr5(13, :)];
ch24_5 = [reconstructed_5(24, :), X_org(24, :), X_observed_snr5(24, :)];
ch13_15 = [reconstructed_15(13, :), X_org(13, :), X_observed_snr15(13, :)];
ch24_15 = [reconstructed_15(24,:), X_org(24,:), X_observed_snr15(24,: )];

% Channel 13 Reconstruction after -5dB noise added
offset1 = max(abs(ch13_5(:)));
disp_eeg(reconstructed_5(13, :),offset1,fs,electrode_name);
title("The Reconstrucrted Signal from -5dB Noise - channel 13");

disp_eeg(X_org(13, :),offset1,fs,electrode_name);
title("The Original Signal - channel 13");

disp_eeg(X_observed_snr5(13, :),offset1,fs,electrode_name);
title("The Noisy Signal -5dB - channel 13");

% Channel13 Reconstruction after -15dB noise added
offset11 = max(abs(ch13_15(:)));
disp_eeg(reconstructed_15(13, :),offset11,fs,electrode_name);
title("The Reconstrucrted Signal from -15dB Noise - channel 13");

disp_eeg(X_org(13, :),offset11,fs,electrode_name);
title("The Original Signal - channel 13");

disp_eeg(X_observed_snr15(13, :),offset11,fs,electrode_name);
title("The Noisy Signal -15dB - channel 13");

%Channel24 Reconstruction after -5dB noise added
offset2 = max(abs(ch24_5(:)));
disp_eeg(reconstructed_5(24,:),offset2,fs,electrode_name);
title("The Reconstrucrted Signal from -5dB noise- channel 24");

disp_eeg(X_org(24,:),offset2,fs,electrode_name);
title("The Original Signal - channel 24");

disp_eeg(X_observed_snr5(24,: ),offset2,fs, electrode_name);
title("The Noisy Signal -5dB - channel 24");


%Channel24 Reconstruction after -15dB noise added
offset21 = max(abs(ch24_15(:)));

disp_eeg(reconstructed_15(24,:),offset21,fs,electrode_name);
title("The Reconstrucrted Signal from -15dB noise - channel 24");

disp_eeg(X_org(24,:),offset21,fs,electrode_name);
title("The Original Signal - channel 24");

disp_eeg(X_observed_snr15(24,: ),offset21,fs, electrode_name);
title("The Noisy Signal -15dB - channel 24");



%% 1.8  RRMSE Calculation

RRMSE = sqrt(sum(sum((X_org-reconstructed_5).^2,1),2)) / sqrt(sum(sum((X_org).^2,1),2));
disp("RRMSE of -5dB SNR")
disp(RRMSE)

RRMSE = sqrt(sum(sum((X_org-reconstructed_15).^2,1),2)) / sqrt(sum(sum((X_org).^2,1),2));
disp("RRMSE of -15dB SNR")
disp(RRMSE)










