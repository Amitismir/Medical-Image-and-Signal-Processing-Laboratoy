clc;
clear;
%%loading data

addpath("E:\6th Semester\MISP Lab\MyLab\LAB2\Lab 2_data\Lab 2_data\Lab2_2")

EEG_1 = load("E:\6th Semester\MISP Lab\MyLab\LAB2\Lab 2_data\Lab 2_data\Lab2_2\NewData1.mat").EEG_Sig;
EEG_2 = load("E:\6th Semester\MISP Lab\MyLab\LAB2\Lab 2_data\Lab 2_data\Lab2_2\NewData4.mat").EEG_Sig;
ElectrodeInfo = load("E:\6th Semester\MISP Lab\MyLab\LAB2\Lab 2_data\Lab 2_data\Lab2_2\Electrodes.mat");


%% Part1
ElecNames=ElectrodeInfo.Electrodes.labels;
fs=250 ;
offset1=max(max(abs(EEG_1)));
offset2=max(max(abs(EEG_2)));
% % ploting eegs 

disp_eeg(EEG_1, offset1,fs,ElecNames);
title("EEG1 Time-domain Representation");
disp_eeg(EEG_2, offset2,fs,ElecNames);
title("EEG4 Time-domain Representation");
xlim([0 5]);

%% Part2
clc;
EEG_1 = EEG_1 - mean(EEG_1,2);
EEG_2 = EEG_2 - mean(EEG_2,2);

[F, W, K] =COM2R(EEG_1, size(EEG_1,1));
[F2, W2, K2] = COM2R(EEG_2, size(EEG_2,1));

EEG_1_ICA_Applied = W * EEG_1;
EEG_2_ICA_Applied = W2 * EEG_2;


%% Part 4 - EEG1 Analysis (21 ICA Components)
% Time domain plot of ICA components
ICA_CompNames = cell(21, 1);
for i = 1:21
    ICA_CompNames{i} = sprintf('IC%d', i);
end

%figure;
DISP_EEG(EEG_1_ICA_Applied(1:21, :), offset1, fs, ICA_CompNames);
title('EEG1 ICA Components (Time Domain)');
xlim([0 5]);

% Frequency domain analysis
figure('Position', [100, 100, 1400, 800]);
for i = 1:21
    [pxx, f] = pwelch(EEG_1_ICA_Applied(i, :), [], [], [], fs);
    subplot(3, 7, i);
    plot(f, 10*log10(pxx));
    title(sprintf('IC%d', i));
    xlabel('Frequency (Hz)');
    xlim([0 80]);
    ylabel('Power (dB/Hz)');
    grid on;
end
sgtitle('EEG1 ICA Components - Frequency Domain');

% Spatial topography (mixing matrix columns)
ElectrodeXs = ElectrodeInfo.Electrodes.X;
ElectrodeYs = ElectrodeInfo.Electrodes.Y;
figure('Position', [100, 100, 1200, 1000]);
for idx = 1:21
    subplot(7, 3, idx);
    plottopomap(ElectrodeXs, ElectrodeYs, ElecNames, W(:, idx));
    title(sprintf('IC%d Topography', idx));
end
sgtitle('EEG1 ICA Components - Spatial Distribution');

%% Part 4 - EEG2 Analysis (21 ICA Components)
% Time domain plot of ICA components
%figure;
DISP_EEG(EEG_2_ICA_Applied(1:21, :), offset2, fs, ICA_CompNames);
title('EEG2 ICA Components (Time Domain)');
xlim([0 5]);

% Frequency domain analysis
figure('Position', [100, 100, 1400, 800]);
for i = 1:21
    [pxx, f] = pwelch(EEG_2_ICA_Applied(i, :), [], [], [], fs);
    subplot(3, 7, i);
    plot(f, 10*log10(pxx));
    title(sprintf('IC%d', i));
    xlabel('Frequency (Hz)');
    xlim([0 80]);
    ylabel('Power (dB/Hz)');
    grid on;
end
sgtitle('EEG2 ICA Components - Frequency Domain');

% Spatial topography (mixing matrix columns)
figure('Position', [100, 100, 1200, 1000]);
for idx = 1:21
    subplot(7, 3, idx);
    plottopomap(ElectrodeXs, ElectrodeYs, ElecNames, W2(:, idx));
    title(sprintf('IC%d Topography', idx));
end
sgtitle('EEG2 ICA Components - Spatial Distribution');
%% PART5
SelSources1 = [2,3,5,7,9,11,13,14,15,17,18,19,20,21];
SelSources2 = [2,6,9,10,14,17,19,20,21];

EEG_1_Denoised = F(:,SelSources1) * EEG_1_ICA_Applied(SelSources1,:);
EEG_2_Denoised = F2(:,SelSources2) * EEG_2_ICA_Applied(SelSources2,:);

%% PART6
ElecNames = ElectrodeInfo.Electrodes.labels; 
offset = max(max(abs(EEG_1_Denoised)));
disp_eeg(EEG_1_Denoised, offset, fs, ElecNames);
xlim([0 5]);
title('EEG1 Signals (Cleaned with ICA) from All Channels');

set(gcf, 'Position', [0, 0, 800, 600]);

ElecNames = ElectrodeInfo.Electrodes.labels; 
offset = max(max(abs(EEG_2_Denoised)));
disp_eeg(EEG_2_Denoised, offset, fs, ElecNames);
xlim([0 5]);
title('EEG4 Signal (Cleaned with ICA) from All Channels');

set(gcf, 'Position', [0, 0, 800, 600]);
%%

function DISP_EEG(Z, offset, feq, ElecName)
     [num_channels, num_samples] = size(Z);
      time = (0:num_samples-1) / feq;
    
    
    figure('Position', [100, 100, 1400, 800]);
    
    % Plot each channel with vertical offset
    hold on;
    
    for i = 1:num_channels
        % Apply offset
        % Subtract mean to center each channel around its offset level
        ch_data = Z(i, :) - mean(Z(i, :)) + (i-1) * offset;
        
        plot(time, ch_data, 'b-', 'LineWidth', 0.5);
        
        % Add channel label on the right side
        text(time(end) + 0.5, (i-1)*offset, ElecName{i}, ...
            'FontSize', 8, 'VerticalAlignment', 'middle');
    end
 
    xlabel('Time (seconds)', 'FontSize', 12);
    ylabel([inputname(1) ' (with offset)'], 'FontSize', 12);
    title('Multichannel EEG Recording', 'FontSize', 14);
    
    % Adjust y-axis limits to show all channels with some padding
    ylim([-offset, (num_channels) * offset]);
    
    % Adjust x-axis limits (add small margin for labels)
    xlim([0, time(end) + 2]);
    
    grid on;
    hold off;
    
    % Optional: Add horizontal lines separating channels
    for i = 0:num_channels
        line([0, time(end)], [i*offset, i*offset], ...
            'Color', [0.8 0.8 0.8], 'LineStyle', ':', 'LineWidth', 0.5);
    end
end 



