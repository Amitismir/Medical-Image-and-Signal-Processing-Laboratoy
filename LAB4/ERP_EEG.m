clear; clc; close all;
% Load the data
data=load("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\LAB4\data\ERP_EEG.mat"); 
ERP_EEG1=data.ERP_EEG;
Fs = 240; % Sampling frequency (Hz)
t = (0:239) / Fs; % Time vector (0 to ~1 second)
total_trials = 2550;

%% Part (a): Average response for N = 100:100:2500

figure;
hold on;
colors = parula(25); % Colormap for 25 lines
c_idx = 1;

for N = 100:100:2500
    % Average across the first N columns (dimension 2)
    avg_response = mean(ERP_EEG1(:, 1:N), 2);

    plot(t, avg_response, 'Color', colors(c_idx,:), 'DisplayName', ['N = ', num2str(N)]);
    c_idx = c_idx + 1;
end
hold off;
title('Part (a): Average ERP for N = 100 to 2500');
xlabel('Time (s)');
ylabel('Amplitude');
% legend('show', 'Location', 'eastoutside'); % Optional: uncomment to see legend
grid on;
%% Part (b): Max absolute amplitude vs N (1 to 2550)
max_amp = zeros(1, total_trials);

for N = 1:total_trials
    % Average the first N trials
    avg_response = mean(ERP_EEG1(:, 1:N), 2);

    % Store the max absolute amplitude of the averaged ERP
    max_amp(N) = max(abs(avg_response));
end

figure;
plot(1:total_trials, max_amp, 'LineWidth', 1.5);
title('Part (b): Maximum Absolute Amplitude vs Number of Trials (N)');
xlabel('Number of Averaged Trials (N)');
ylabel('Max Absolute Amplitude');
grid on;


%% Part (c): Plot RMSE between i-th and (i-1)-th average templates

N_total_in_file = 2550; 
n_values_for_rmse = 100:100:2500;
rmse_vals = zeros(1, length(n_values_for_rmse));

for i = 2:length(n_values_for_rmse)
    n = n_values_for_rmse(i);       % Current number of trials for averaging
    n_minus_1 = n_values_for_rmse(i-1); % Previous number of trials for averaging

    current_n_trials_data = ERP_EEG1(:, 1:n);
    previous_n_trials_data = ERP_EEG1(:, 1:n_minus_1);
    avg_erpn = mean(current_n_trials_data, 2);
    avg_erpn_minus_1 = mean(previous_n_trials_data, 2);

    % Calculate Root Mean Square Error (RMSE) between the two average templates.
    rmse_vals(i) = sqrt(mean((avg_erpn - avg_erpn_minus_1).^2));
end

figure;
plot(n_values_for_rmse(2:end), rmse_vals(2:end), '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Number of Trials Averaged (N)');
ylabel('RMSE between N and (N-1) Averages');
title('RMSE between Consecutive Average Templates (N vs N-1)');
grid on;

%% part(d)
%N0​≈500: This is often considered the "knee” of the curve. It is a good, efficient choice because you achieve a reasonably 
% stable template with a moderate number of trials. This is often “better” if you want to minimize recording time for the subject.
%N0​≈1000: This is a more conservative choice. While it takes longer to record, the RMSE is lower and the curve is flatter, meaning the
% resulting P300 template is more robust and less susceptible to minor noise variations.

%% part(e)
N_0 = 500;  
N_2550 = 2550;
N_0_3 = 150;

avg_erp_N_0 = mean(ERP_EEG1(:, 1:N_0), 2); 
avg_erp_N_0_3 = mean(ERP_EEG1(:, 1:N_0_3), 2); 
avg_erp_2550 = mean(ERP_EEG1, 2); 

random_N_0_trials = ERP_EEG1(:, randperm(2550, N_0));
avg_erp_random_N_0_trials = mean(random_N_0_trials, 2);  

random_N_0_3_trials = ERP_EEG1(:, randperm(2550, N_0_3));
avg_erp_random_N_0_3_trials = mean(random_N_0_3_trials, 2);  

figure;
hold on;
plot(t, avg_erp_N_0, 'LineWidth', 2, 'DisplayName', 'Average for N = N0 = 500');
plot(t, avg_erp_2550, 'LineWidth', 2, 'DisplayName', 'Average for N = 2550');
plot(t, avg_erp_N_0_3, 'LineWidth', 2, 'DisplayName', 'Average for N = N0/3 = 140');

plot(t, avg_erp_random_N_0_trials, 'LineWidth', 2, 'DisplayName', 'Average for Random N = random N0= 500');
plot(t, avg_erp_random_N_0_3_trials, 'LineWidth', 2, 'DisplayName', 'Average for Random N = random N0/3 = 150');


ylabel('Amplitude (\muV)');
title('Comparison of Different Averaged Responses');
legend;
grid on;