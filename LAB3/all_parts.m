%%PART1
% time domain_ plot 
mecg1=load("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\data\mecg1.dat");
fecg1=load("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\data\fecg1.dat");
noise1=load("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\data\noise1.dat");

ECG_signal=mecg1+fecg1+noise1;

Fs=256;
t=(0:length(ECG_signal)-1/Fs);
figure
plot(t, ECG_signal, 'r', 'LineWidth', 1.2);
xlabel('Time (seconds)');
ylabel('Amplitude (mV)');
xlim([0, max(t)]);
title('Combined ECG Signal (Time Domain)');
grid on;
%%
figure;
subplot(3,1,1);
plot(t, mecg1); title('Maternal ECG'); ylabel('mV');
xlim([0, max(t)]);

subplot(3,1,2);
plot(t, fecg1); title('Fetal ECG'); ylabel('mV');
xlim([0, max(t)]);

subplot(3,1,3);
plot(t, noise1); title('Noise'); ylabel('mV'); xlabel('Time (s)');
xlim([0, max(t)]);
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
X_data=load("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\data\X.dat");
addpath("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\matlab");
plot3ch(X_data)

%% P2 3D Visualization of Singular Vectors with Corresponding Singular Values

[U, S, V] = svd(X_data, 'econ');

for v = 1:size(V, 2)
    plot3dv(V(:, v), S(v, v));  % Custom function provided in lab
    hold on;
end

title('3D Visualization of Singular Vectors with Corresponding Singular Values');
xlabel('X Axis');
ylabel('Y Axis');
zlabel('Singular Value (S)');
grid on;


savefig('svd_singular_vectors.fig');

% Save U, S, V to MAT file
save('svd_decomposition_results.mat', 'U', 'S', 'V');

%%  p3


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
    xlim([0, max(t)]);
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
    plot(X_reconstructed(:, i), 'LineWidth', 0.5);
    title(['Reconstructed Channel (Column ', num2str(i), ' of U)']);
    xlabel('Time Samples');
    ylabel('Amplitude');
    xlim([0, max(t)]);
    grid on;
end

%% Part 3.1
addpath("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\data")
addpath("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\matlab")
X_mat = load("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\data\X.dat");
X_mat_t = X_mat';
[W_mat, Z_mat] = ica(X_mat_t);
A = inv(W_mat);
save('Z_mat.mat','Z_mat');
save('W_mat.mat','W_mat');
save('A_mat.mat','A');
%% Part 3.2.1
clc;
plot3ch(X_mat);
%% Part 3.2.2
clc;
plot3dv(A(:,1));
plot3dv(A(:,2));
plot3dv(A(:,3));
savefig(figure(1),'part3.fig');
%% Part 3.3.1
clc;
figure;
for i = 1:3
    subplot(3,1,i);
    plot( Z_mat(i,:));
    title(['Channel', num2str(i)]);
    xlabel('Time (samples)');
    ylabel('Amplitude(a.u.)');
    xlim([0 2560]);
    grid on;
end
%% Part 3.3.2
clc;
A_modified = A;
A_modified(:,1:2) = 0;
X_recons = A_modified * Z_mat;
%% Part 3.4
clc;
plot3ch(X_recons', 256, "Reconstructed ECG")
%% Part 4.1
figure;
plot3(X_mat(:,1), X_mat(:,2), X_mat(:,3), '.m');
hold all;
plot3(X_reconstructed(:, 1), X_reconstructed(:, 2), X_reconstructed(:, 3), '+'); % This is from SVD
hold all;
plot3(X_recons(1,: ), X_recons(2, :), X_recons(3, :), '*'); % This is from ICA
xlabel('Ch1');
ylabel('Ch2');
zlabel('Ch3');
title('Comparison of scatter plots');
grid on;
plot3dv(A(:,1), 100, 'black');
plot3dv(A(:,2), 100, 'black');
plot3dv(A(:,3), 100, 'black');
plot3dv(V(:,1), 100, 'blue'); 
plot3dv(V(:,2), 100, 'blue'); 
plot3dv(V(:,3), 100, 'blue');

norm_A = [norm(A(:,1)),norm(A(:,2)),norm(A(:,3))];
angles_A = [dot(A(:,1),A(:,2))/norm_A(1)/norm_A(2),dot(A(:,1),A(:,3))/norm_A(1)/norm_A(3),dot(A(:,2),A(:,3))/norm_A(2)/norm_A(3)];
angles_A = acosd(angles_A);

norm_V = [norm(V(:,1)),norm(V(:,2)),norm(V(:,3))];
angles_V = [dot(V(:,1),V(:,2))/norm_V(1)/norm_V(2),dot(V(:,1),V(:,3))/norm_V(1)/norm_V(3),dot(V(:,2),V(:,3))/norm_V(2)/norm_V(3)];
angles_V = acosd(angles_V);
%% Part 4.2
clc;
FECG2 = load("E:\6th Semester\MISP Lab\MyLab\LAB3\Lab 3\data\fecg2.dat");
figure;
subplot(3,1,1);
plot(t, FECG2);
hold on;
title('Original ECG Fetus Signal ');
xlabel("Time(us)");
xlim([0 max(t)]);
grid on;

subplot(3,1,2);
plot(t, X_reconstructed(:,1));
hold on;
title('The Reconstructed Signal Using SVD');
xlabel("Time(us)");
xlim([0 max(t)]);
grid on;

subplot(3,1,3);
plot(t, X_recons(1, :));
hold on;
title('The Reconstructed Signal Using ICA');
xlabel("Time(us)");
xlim([0 max(t)]);
grid on;
%% 4.4
clc;
corr_svd = corrcoef(X_reconstructed(:,1), FECG2(:,1));
corr_ica = corrcoef(X_recons(1,:)', FECG2(:,1));




