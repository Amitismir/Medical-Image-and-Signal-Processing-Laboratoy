%%PART1
% time domain_ plot 
mecg1=load("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\LAB3\Data\mecg1.dat");
fecg1=load("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\LAB3\Data\fecg1.dat");
noise1=load("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\LAB3\Data\noise1.dat");

ECG_signal=mecg1+fecg1+noise1;

Fs=256;
t=(0:length(ECG_signal)-1/Fs);
figure
plot(t, ECG_signal, 'r', 'LineWidth', 1.2);
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
title('Combined ECG Signal (Time Domain)');
grid on;
%%
figure;
subplot(3,1,1);
plot(t, mecg1); title('Maternal ECG'); ylabel('mV');

subplot(3,1,2);
plot(t, fecg1); title('Fetal ECG'); ylabel('mV');

subplot(3,1,3);
plot(t, noise1); title('Noise'); ylabel('mV'); xlabel('Time (s)');
%% PART 2
fs = 256;

[pxx_m, f_m] = pwelch(mecg1, [], [], [], fs);
[pxx_f, f_f] = pwelch(fecg1, [], [], [], fs);
[pxx_n, f_n] = pwelch(noise1, [], [], [], fs);

figure;

subplot(3,1,1)
plot(f_m,10*log10(pxx_m))
title('Maternal ECG Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Power (dB/Hz)')
xlim([0 80])
grid on

subplot(3,1,2)
plot(f_f,10*log10(pxx_f))
title('Fetal ECG Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Power (dB/Hz)')
xlim([0 80])
grid on

subplot(3,1,3)
plot(f_n,10*log10(pxx_n))
title('Noise Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Power (dB/Hz)')
xlim([0 80])
grid on

%% p3_ Mean and Variance

mean_m = mean(mecg1);
var_m  = var(mecg1);

mean_f = mean(fecg1);
var_f  = var(fecg1);

mean_n = mean(noise1);
var_n  = var(noise1);

fprintf("Maternal ECG:   Mean = %.4f   | Variance = %.4f\n", mean_m, var_m);
fprintf("Fetal ECG:      Mean = %.4f   | Variance = %.4f\n", mean_f, var_f);
fprintf("Noise:          Mean = %.4f   | Variance = %.4f\n", mean_n, var_n);

figure;
% Maternal ECG
subplot(3,1,1);
histogram(mecg1, 100, 'Normalization','pdf')
title('Histogram and PDF of Maternal ECG')
xlabel('Amplitude (mV)')
ylabel('PDF')
grid on

% Fetal ECG
subplot(3,1,2);
histogram(fecg1, 100, 'Normalization','pdf')
title('Histogram and PDF of Fetal ECG')
xlabel('Amplitude (mV)')
ylabel('PDF')
grid on

% Noise
subplot(3,1,3);
histogram(noise1, 100, 'Normalization','pdf')
title('Histogram and PDF of Noise')
xlabel('Amplitude (mV)')
ylabel('PDF')
grid on
k_m = kurtosis(mecg1);
k_f = kurtosis(fecg1);
k_n = kurtosis(noise1);

fprintf('Kurtosis (Maternal ECG): %.4f\n', k_m);
fprintf('Kurtosis (Fetal ECG): %.4f\n', k_f);
fprintf('Kurtosis (Noise): %.4f\n', k_n);

%% SVD  PART2
%%p1
X_data=load("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\LAB3\Data\X.dat");
addpath("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\LAB3\Data\plot3ch.m");
figure;
plot3ch(X_data,fs,"X_data")

%% P2 3D Visualization of Singular Vectors with Corresponding Singular Values

figure;

for v = 1:size(V, 2)
    plot3dv(V(:, v), S(v, v));  % Custom function provided in lab
    hold on;
end

title('3D Visualization of Singular Vectors with Corresponding Singular Values');
xlabel('X Axis');
ylabel('Y Axis');
zlabel('Singular Value (S)');
grid on;

 Save figure
savefig('svd_singular_vectors.fig');

% Save U, S, V to MAT file
save('svd_decomposition_results.mat', 'U', 'S', 'V');

%%  p3
[U, S, V] = svd(X_data, 'econ');

% 3. Extract Singular Values (Diagonal elements of S)
singular_values = diag(S);

% 4. Plotting the Eigenspectrum
figure('Name', 'SVD Analysis - Eigenspectrum', 'Color', 'w');
stem(singular_values, 'filled', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on;
title('Eigenspectrum (Singular Values of S)');
xlabel('Component Index');
ylabel('Singular Value Magnitude');
set(gca, 'FontSize', 12);

% 5. Plotting the Principal Components (Columns of U)
figure('Name', 'Principal Components (U Matrix)', 'Color', 'w');

for i = 1:3
    subplot(3, 1, i);
    plot(U(:, i), 'LineWidth', 1);
    title(['Principal Component (Column ', num2str(i), ' of U)']);
    xlabel('Time Samples');
    ylabel('Amplitude');
    grid on;
end
fprintf('Size of U: %d x %d\n', size(U,1), size(U,2));
fprintf('Size of S: %d x %d\n', size(S,1), size(S,2));
fprintf('Size of V: %d x %d\n', size(V,1), size(V,2));
%%
%% Intuition and Conceptual Understanding of SVD
% Why do we use SVD?
% Imagine the recorded ECG signal is a mixture of several sources:
% maternal ECG, fetal ECG, and some background noise.
% SVD helps us "unmix" this composite signal and reveal the underlying components.

% (a) Matrix S – Energy and Importance
% The diagonal elements of matrix S are the singular values.
% Each singular value represents the amount of energy, or variance,
% carried by its corresponding component.
% In ECG signals, the first singular value is usually much larger
% because the maternal heartbeat dominates the overall energy.
% Subsequent singular values correspond to the fetal heartbeat and then noise.
% If you plot the singular values (the Eigenspectrum),
% you will typically see a sharp drop after the first few values.
% This drop indicates that the main signal information
% is concentrated in the first few components,
% while smaller singular values mostly represent noise.

% (b) Matrix V – Feature Space
% The columns of V are orthogonal vectors that describe directions
% of maximal variance in the data.
% These are equivalent to eigenvectors of the covariance matrix (C = Y' * Y).
% Intuitively, V tells us how correlated or linearly dependent the signals are
% with respect to each other—it describes the "feature space" structure.

% (c) Matrix U – Independent Waveforms
% The columns of U are new versions of the input signals,
% expressed in an orthogonal basis where each component is independent.
% In practice, when you perform SVD on the composite ECG signal,
% the first column of U often represents a cleaned version
% of the maternal heartbeat waveform,
% while subsequent columns may represent the fetal ECG or noise.
% Therefore, SVD is a useful technique for isolating the main physiological
% components and removing interference or noise
% by reconstructing the signal using only the strongest singular values.
%% p4
 S_new = S;
 S_new(:, 1) = 0;
 S_new(:, 3) = 0;
 X_reconstructed = U * S_new* V';
 figure('Name', 'Reconstructed Channel', 'Color', 'w');

for i = 1:3
    subplot(3, 1, i);
    plot(U(:, i), 'LineWidth', 0.5);
    title(['Reconstructed Channel (Column ', num2str(i), ' of U)']);
    xlabel('Time Samples');
    ylabel('Amplitude');
    grid on;
end

