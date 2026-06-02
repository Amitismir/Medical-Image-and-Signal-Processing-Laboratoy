%% 2.1
clc;
load("E:\6th Semester\MISP Lab\MyLab\LAB5\Lab 5_data\Lab 5_data\n_422.mat");
fs = 250;
t = (0:length(n_422)-1)/fs;


normal_seg = n_422(1:1+fs*10);
abnormal_seg = n_422(11442 : 11442 + fs*10);
[b,a] = butter(4,[0.5 45]/ (fs/2),"bandpass");

Filtered_normal = filtfilt(b,a,normal_seg);
Filtered_abnormal = filtfilt(b,a,abnormal_seg);

t = linspace(0, length(normal_seg)/fs, length(normal_seg)); 

[p_normal, f_normal] = pwelch(Filtered_normal, [], [], [], fs);
[p_abnormal, f_abnormal] = pwelch(Filtered_abnormal, [], [], [], fs);


%% 2.3
fs = 250;
window_length = 10 * fs;
overlap = 5 * fs;
step = window_length - overlap;

num_windows = floor((length(n_422) - window_length)/step) + 1;

event_labels = [
    1       10711    1;   % Normal
    10711   11211    3;   % VT
    11211   11442    1;   % Normal
    11442   59711    3;   % VT
    59711   61288    4;   % Noise
    61288   75000    2    % VFIB
];

true_labels = zeros(num_windows, 1);
end_samples = zeros(num_windows, 1);

for i = 1:num_windows
    win_start = (i-1)*step + 1;
    win_end = win_start + window_length - 1;
    end_samples(i) = win_end;
    
    label = 0;  
    for j = 1:size(event_labels, 1)
        range_start = event_labels(j, 1);
        range_end = event_labels(j, 2);
        range_label = event_labels(j, 3);
        
        if win_start >= range_start && win_end <= range_end
            label = range_label;
            break;
        elseif (win_start < range_end && win_end > range_start)
            label = 0;
        end
    end
    
    true_labels(i) = label;
end




%% 2.4 and 2.5
features = zeros(length(true_labels), 5);
num_windows = length(true_labels);

for i = 1:num_windows
    start_idx = (i - 1) * overlap + 1; 
    end_idx = start_idx + window_length - 1;

    if end_idx > length(n_422)
        break;
    end

    window = n_422(start_idx:end_idx);
    features(i, :) = extract_frequency_features(window, fs);
end

    features_healthy = features(true_labels == 1, :);
    features_vfib = features(true_labels==2, :);

    feature_names = {'Dominant Freq', 'Mean Freq', 'Median Freq', 'Spectral Entropy', 'Band Power'};



%% 2.6 Alram Vector
valid_indices = (true_labels==1) | (true_labels == 2);

GT_h = true_labels;
GT_h(~valid_indices) = 0;
[alarm_n422_h,t] = va_detect_424(n_422,fs, valid_indices);

valid_indices_2 = [find(true_labels==1) ; find(true_labels == 2)];
GT=GT_h(valid_indices_2);
alarm_n422=alarm_n422_h(valid_indices_2);
%% 2.7 
ground_truth = GT;
predicted = alarm_n422;

figure;
confusionchart(ground_truth, predicted, ...
    'RowSummary','row-normalized', ...
    'ColumnSummary','column-normalized', ...
    'Title','Confusion Matrix');

confMat = confusionmat(ground_truth, predicted);
disp('Confusion Matrix:');
disp(confMat);

accuracy = sum(diag(confMat)) / sum(confMat(:));
fprintf('Overall Accuracy: %.2f%%\n\n', accuracy * 100);

class_labels = [0 1 2];
n_classes = numel(class_labels);

for i = 1:n_classes
    class = class_labels(i);

    TP = confMat(class+1, class+1);  
    FN = sum(confMat(class+1, :)) - TP;  
    FP = sum(confMat(:, class+1)) - TP; 
    TN = sum(confMat(:)) - TP - FN - FP;

    sensitivity = TP / (TP + FN);
    specificity = TN / (TN + FP);

    fprintf('Class %d:\n', class);
    fprintf('  Sensitivity (Recall): %.2f%%\n', sensitivity * 100);
    fprintf('  Specificity: %.2f%%\n\n', specificity * 100);
end


%% Functions
function feat = extract_frequency_features(signal, fs)
    [pxx, f] = pwelch(signal, [], [], [], fs);
    pxx_norm = pxx / sum(pxx + eps);
      
    [~, idx] = max(pxx);
    DF = f(idx); 
    MF = sum(f .* pxx) / sum(pxx);
    cumulative_power = cumsum(pxx);
    total_power = cumulative_power(end);
    idx_median = find(cumulative_power >= total_power/2, 1);
    MedF = f(idx_median);
    SpecEnt = -sum(pxx_norm .* log2(pxx_norm + eps));
    band_idx = f >= 0.5 & f <= 40;
    BandPower = sum(pxx(band_idx));
    feat = [DF, MF, MedF, SpecEnt, BandPower];
end

function [alarm,t] = va_detect(ecg_data,Fs, GT)

frame_sec = 10;  
overlap = 0.5;    

ecg_data = ecg_data(:);  

frame_length = round(frame_sec*Fs);  % length of each data frame (samples)
frame_step = round(frame_length*(1-overlap));  % amount to advance for next data frame
ecg_length = length(ecg_data);  % length of input vector
frame_N = floor((ecg_length-(frame_length-frame_step))/frame_step); % total number of frames
alarm = zeros(frame_N,1);	% initialize output signal to all zeros
t = ([0:frame_N-1]*frame_step+frame_length)/Fs;
DF_h = [];
MedF_h = [];
for i = 1:frame_N
    seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));
    if (GT(i) == 1)
        [pxx, f] = pwelch(seg, [], [], [], Fs);    
        [~, idx] = max(pxx);
        DF = f(idx);
        DF_h = [DF_h, DF];
    
        cumulative_power = cumsum(pxx);
        total_power = cumulative_power(end);
        idx_median = find(cumulative_power >= total_power/2, 1);
        MedF = f(idx_median);
        MedF_h = [MedF_h, MedF];
    end
end

AV_DF=mean(DF_h);
AV_medF=mean(MedF_h);

for i = 1:frame_N
    if (GT(i) == 1)
        seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));

        [pxx, f] = pwelch(seg, [], [], [], Fs);    
        [~, idx] = max(pxx);
        DF = f(idx);
        cumulative_power = cumsum(pxx);
        total_power = cumulative_power(end);
        idx_median = find(cumulative_power >= total_power/2, 1);
        MedF = f(idx_median);

        if DF>AV_DF && MedF>AV_medF
            alarm(i) = 2;
        elseif DF<AV_DF && MedF<AV_medF
            alarm(i) = 1;
        else 
            alarm(i) = 0;
        end
    else
        alarm(i) = 0;
    end
end
end


function featArray = extract_morphological_features(signal)
 
    signal = signal(:);

    maxAmp = max(signal);
 
    minAmp = min(signal);
    
    ptpAmp = maxAmp - minAmp;
    
    zeroCrossings = sum(diff(sign(signal)) ~= 0);
    
    ampVariance = var(signal);
    [pks, ~] = findpeaks(signal, 'MinPeakHeight', 0.5 * maxAmp, 'MinPeakDistance', round(0.25 * 250)); 
    if isempty(pks)
        avgRPeakAmp = NaN; 
    else
        avgRPeakAmp = mean(pks);
    end

    featArray = [maxAmp, minAmp, ptpAmp, zeroCrossings, ampVariance, avgRPeakAmp];
end


function [alarm,t] = va_detect_2(ecg_data,Fs, GT)

frame_sec = 10;  
overlap = 0.5;    

ecg_data = ecg_data(:);  

frame_length = round(frame_sec*Fs);  % length of each data frame (samples)
frame_step = round(frame_length*(1-overlap));  % amount to advance for next data frame
ecg_length = length(ecg_data);  % length of input vector
frame_N = floor((ecg_length-(frame_length-frame_step))/frame_step); % total number of frames
alarm = zeros(frame_N,1);	% initialize output signal to all zeros
t = ([0:frame_N-1]*frame_step+frame_length)/Fs;
ZC_h = [];
RPeak_h = [];
for i = 1:frame_N
    seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));
    if (GT(i) == 1)
        ZC = sum(diff(sign(seg)) ~= 0);
        ZC_h = [ZC_h, ZC];
        
        [pks, ~] = findpeaks(seg, 'MinPeakHeight', 0.5 * (max(seg)), 'MinPeakDistance', round(0.25 * 250)); 
        if isempty(pks)
            RPeak = NaN;  % In case no peaks are found
        else
            RPeak = mean(pks);
        end
        RPeak_h = [RPeak_h, RPeak];
    end
end

AV_ZC = mean(ZC_h);
AV_RPeak = mean(RPeak_h);

for i = 1:frame_N
    if (GT(i) == 1)
        seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));

        ZC = sum(diff(sign(seg)) ~= 0);
        
        [pks, ~] = findpeaks(seg, 'MinPeakHeight', 0.5 * (max(seg)), 'MinPeakDistance', round(0.25 * 250)); 
        if isempty(pks)
            RPeak = NaN;  % In case no peaks are found
        else
            RPeak = mean(pks);
        end

        if ZC > AV_ZC && RPeak < AV_RPeak
            alarm(i) = 2;
        elseif ZC < AV_ZC && RPeak > AV_RPeak
            alarm(i) = 1;
        else 
            alarm(i) = 0;
        end
    else
        alarm(i) = 0;
    end
end
end


function [alarm,t] = va_detect_424(ecg_data,Fs, GT)

frame_sec = 10;  
overlap = 0.5;    

ecg_data = ecg_data(:);  

frame_length = round(frame_sec*Fs);  % length of each data frame (samples)
frame_step = round(frame_length*(1-overlap));  % amount to advance for next data frame
ecg_length = length(ecg_data);  % length of input vector
frame_N = floor((ecg_length-(frame_length-frame_step))/frame_step); % total number of frames
alarm = zeros(frame_N,1);	% initialize output signal to all zeros
t = ([0:frame_N-1]*frame_step+frame_length)/Fs;

MF_h = [];
SE_h = [];
for i = 1:frame_N
    seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));
    [pxx, f] = pwelch(seg, [], [], [], 250);
    pxx_norm = pxx / sum(pxx + eps);

    if (GT(i) == 1)
        MF = sum(f .* pxx) / sum(pxx);
        MF_h = [MF_h, MF];
    
        SE = -sum(pxx_norm .* log2(pxx_norm + eps));
        SE_h = [SE_h, SE];
    end
end

AV_MF = mean(MF_h);
AV_SE = mean(SE_h);

for i = 1:frame_N
    if (GT(i) == 1)
        seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));
        [pxx, f] = pwelch(seg, [], [], [], 250);
        pxx_norm = pxx / sum(pxx + eps);

        MF = sum(f .* pxx) / sum(pxx);
        SE = -sum(pxx_norm .* log2(pxx_norm + eps));

        if MF > AV_MF && SE > AV_SE
            alarm(i) = 2;
        elseif MF < AV_MF && SE < AV_SE
            alarm(i) = 1;
        else 
            alarm(i) = 0;
        end
    else
        alarm(i) = 0;
    end
end
end

function [alarm,t] = va_detect_2_424(ecg_data,Fs, GT)

frame_sec = 10;  
overlap = 0.5;    

ecg_data = ecg_data(:);  

frame_length = round(frame_sec*Fs);  % length of each data frame (samples)
frame_step = round(frame_length*(1-overlap));  % amount to advance for next data frame
ecg_length = length(ecg_data);  % length of input vector
frame_N = floor((ecg_length-(frame_length-frame_step))/frame_step); % total number of frames
alarm = zeros(frame_N,1);	% initialize output signal to all zeros
t = ([0:frame_N-1]*frame_step+frame_length)/Fs;
ZC_h = [];
AV_h = [];

for i = 1:frame_N
    seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));

    if (GT(i) == 1)
        ZC = sum(diff(sign(seg)) ~= 0);
        ZC_h = [ZC_h, ZC];
        
        AV = var(seg);
        AV_h = [AV_h, AV];
    end
end

AV_ZC = mean(ZC_h);
AV_AV = mean(AV_h);

for i = 1:frame_N
    if (GT(i) == 1)
        seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));

        ZC = sum(diff(sign(seg)) ~= 0);
        
        AV = var(seg);

        if ZC > AV_ZC && AV < AV_AV
            alarm(i) = 2;
        elseif ZC < AV_ZC && AV > AV_AV
            alarm(i) = 1;
        else 
            alarm(i) = 0;
        end
    else
        alarm(i) = 0;
    end
end
end

