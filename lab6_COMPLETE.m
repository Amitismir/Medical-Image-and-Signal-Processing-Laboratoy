
clc; clear; close all;
addpath("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\MISP Lab - Source\Lab 6\Codes and Data\");

%% Part (a): Model Configuration and Setup 
load("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\MISP Lab - Source\Lab 6\Codes and Data\ElecPosXYZ.mat");

ModelParams.R      = [9.2, 8.5, 8.0];      
ModelParams.Sigma  = [3.3e-3, 8.25e-5, 3.3e-3]; 
ModelParams.Lambda = [0.5979, 0.2037, 0.0237];
ModelParams.Mu     = [0.6342, 0.9364, 1.0362];

gridRes = 1;   

[LocMat, GainMat] = ForwardModel_3shell(gridRes, ModelParams);

nElec    = 21;          
nSources = size(LocMat, 2);  
%% Part (B): Visualize Source Grid and Electrode Positions
fig_main = figure('Name', 'Source Grid and Electrode Positions');
scatter3(LocMat(1,:), LocMat(2,:), LocMat(3,:), 5, [0.2 0.6 1], 'filled');
hold on;

ElecXYZ = zeros(nElec, 3);
ElecNames = cell(nElec, 1);
for k = 1:nElec
    ElecXYZ(k,:) = ElecPos{k}.XYZ * ModelParams.R(3);
    ElecNames{k} = ElecPos{k}.Name;
    text(ElecXYZ(k,1), ElecXYZ(k,2), ElecXYZ(k,3), ElecNames{k}, 'FontSize', 8, 'FontWeight', 'bold');
end
scatter3(ElecXYZ(:,1), ElecXYZ(:,2), ElecXYZ(:,3), 80, 'r', 'filled');

xlabel('X (cm)'); ylabel('Y (cm)'); zlabel('Z (cm)');
title('Dipole Source Grid and Electrode Positions on Scalp');
legend({'Source Grid', 'Electrodes'}, 'Location', 'best');
grid on; axis equal;
view(45, 30);

%%  Part (پ): Select Dipole Based on Location Criteria (3 Cases) -
dipoleCases = {'deep', 'cortex_central', 'cortex_temporal'};
caseNames = {'Deep Dipole (بخش عمق)', 'Central Cortex Dipole (بخش مرکزی)', 'Temporal Cortex Dipole (بخش تمپورال)'};

allResults = struct();

for caseIdx = 1:length(dipoleCases)
    dipoleType = dipoleCases{caseIdx};

   % Select dipole matching the location criterion
    foundDipole = false;
    attempts = 0;
    while ~foundDipole && attempts < 1000
        rndIdx = randi(nSources);
        dipNorm = norm(LocMat(:, rndIdx));
        r1 = ModelParams.R(1);
        r3 = ModelParams.R(3);
        
        switch dipoleType
            case 'cortex_central'  % Central cortex (near vertex)
                foundDipole = (dipNorm < r1) && (dipNorm > 0.8*r1) && ...
                              (abs(LocMat(1,rndIdx)) < 3) && (abs(LocMat(2,rndIdx)) < 3);
                              
            case 'cortex_temporal'  % Temporal cortex (lateral)
                foundDipole = (dipNorm < r1) && (dipNorm > 0.8*r1) && ...
                              (abs(LocMat(2,rndIdx)) > 0.7*r3) && ...
                              (abs(LocMat(1,rndIdx)) < 4);
                              
            case 'deep'  % Deep dipole (subcortical)
                foundDipole = dipNorm < 0.2 * r1;
        end
        attempts = attempts + 1;
    end
    
    if ~foundDipole
        warning('Could not find dipole for case: %s', dipoleType);
        continue;
    end
    
    % Dipole orientation = radial unit vector
    dipOrient = LocMat(:, rndIdx) / dipNorm;
    
    fprintf('Selected dipole index: %d\n', rndIdx);
    fprintf('Dipole location (x,y,z): [%.2f, %.2f, %.2f] cm\n', LocMat(1,rndIdx), LocMat(2,rndIdx), LocMat(3,rndIdx));
    fprintf('Dipole depth (radius): %.4f cm\n', dipNorm);
    fprintf('Dipole orientation (radial): [%.3f, %.3f, %.3f]\n', dipOrient(1), dipOrient(2), dipOrient(3));
    
    %% Part (ت): Simulate EEG Measurements
    load("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\MISP Lab - Source\Lab 6\Codes and Data\Interictal.mat");
    
    % Lead field for selected dipole
    dipCols = (rndIdx-1)*3 + (1:3);
    GainDip = GainMat(:, dipCols);
    
    % Simulate EEG potentials
    ElecPot = GainDip * (Interictal(1,:) .* dipOrient);
    
    % Display EEG signals
    figure('Name', sprintf('EEG Signals - %s', caseNames{caseIdx}));
    addpath("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\MISP Lab - Source\Lab 2\Codes\");
    disp_eeg(ElecPot, 15, 256, ElecNames);
    title(sprintf('Simulated EEG Potentials - %s', caseNames{caseIdx}));
    
    %% Part (ث): Peak-Epoch Averaging and Scalp Potential Map
    figure('Name', sprintf('Mean Scalp Potential - %s', caseNames{caseIdx}));
    
    % Compute average potential around peaks
    peakMean = zeros(nElec, 1);
    for ch = 1:nElec
        sig = ElecPot(ch,:);
        [pks, pkLocs] = findpeaks(sig, 'MinPeakProminence', 0.9*max(sig));
        
        if ~isempty(pkLocs)
            halfWin = 3;
            validPks = pkLocs(pkLocs > halfWin & pkLocs <= length(sig)-halfWin);
            if ~isempty(validPks)
                epochs = zeros(length(validPks), 2*halfWin+1);
                for p = 1:length(validPks)
                    epochs(p,:) = sig(validPks(p)-halfWin : validPks(p)+halfWin);
                end
                peakMean(ch) = mean(epochs, 'all');
            end
        end
    end
    
    addpath("C:\Users\rakyn\OneDrive\Desktop\TERM6\MISP_LAB\LAB6\Codes and Data\");
    Display_Potential_3D(ModelParams.R(3), peakMean);
    title(sprintf('Scalp Potential Map - %s', caseNames{caseIdx}));
    
    %% Part (ج): Minimum Norm Estimation (MNE) Source Localization 
    alpha = 0.1;  
    G = GainMat;
    M = ElecPot;
    
    
    Q_hat = G' * ((G*G' + alpha*eye(nElec)) \ M);
 
    qNorms = zeros(nSources, 1);
    for s = 1:nSources
        idx3 = (s-1)*3 + (1:3);
        qNorms(s) = norm(Q_hat(idx3, :));
    end
    
    [~, estIdx] = max(qNorms);
    estLoc = LocMat(:, estIdx);
    estNorm = norm(estLoc);
    estOrient = estLoc / estNorm;
   
    locError = norm(estLoc - LocMat(:, rndIdx));
    angleError = acosd(abs(dot(estOrient, dipOrient)));
    
    fprintf('\n--- MNE Localization Results ---\n');
    fprintf('True dipole location: [%.2f, %.2f, %.2f] cm\n', LocMat(1,rndIdx), LocMat(2,rndIdx), LocMat(3,rndIdx));
    fprintf('Estimated location:   [%.2f, %.2f, %.2f] cm\n', estLoc(1), estLoc(2), estLoc(3));
    fprintf('Location error: %.4f cm\n', locError);
    fprintf('Orientation error: %.2f degrees\n', angleError);
    
    %%  Display Real vs Estimated Dipole Vectors
    fig_compare = figure('Name', sprintf('Real vs Estimated - %s', caseNames{caseIdx}));
    scatter3(LocMat(1,:), LocMat(2,:), LocMat(3,:), 5, [0.2 0.6 1], 'filled');
    hold on;
    
    % Plot electrodes
    for k = 1:nElec
        scatter3(ElecXYZ(k,1), ElecXYZ(k,2), ElecXYZ(k,3), 60, 'r', 'filled');
        text(ElecXYZ(k,1), ElecXYZ(k,2), ElecXYZ(k,3), ElecNames{k}, 'FontSize', 8);
    end
    
    % Plot true dipole (green)
    quiver3(LocMat(1,rndIdx), LocMat(2,rndIdx), LocMat(3,rndIdx), ...
            dipOrient(1), dipOrient(2), dipOrient(3), ...
            'g', 'LineWidth', 3, 'MaxHeadSize', 1.5);
    scatter3(LocMat(1,rndIdx), LocMat(2,rndIdx), LocMat(3,rndIdx), 100, 'g', 'filled');
    
    % Plot estimated dipole (magenta)
    quiver3(estLoc(1), estLoc(2), estLoc(3), ...
            estOrient(1), estOrient(2), estOrient(3), ...
            'm', 'LineWidth', 3, 'MaxHeadSize', 1.5);
    scatter3(estLoc(1), estLoc(2), estLoc(3), 100, 'm', 'filled');
    
    % Plot origin
    scatter3(0, 0, 0, 100, 'r', '*', 'LineWidth', 2);
    
    % Add error information text
    textStr = sprintf('Location Error: %.4f cm\nOrientation Error: %.2f°', locError, angleError);
    text(-12, -12, 12, textStr, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 10, 'FontWeight', 'bold');
    
    xlabel('X (cm)'); ylabel('Y (cm)'); zlabel('Z (cm)');
    title(sprintf('Real vs Estimated Dipole Vectors - %s', caseNames{caseIdx}));
    legend({'Source Grid', 'Electrodes', 'True Dipole', '', 'Estimated Dipole', '', 'Origin'}, 'Location', 'best');
    grid on; axis equal;
    view(45, 30);
    
    %% Store results
    allResults(caseIdx).case = caseNames{caseIdx};
    allResults(caseIdx).trueLoc = LocMat(:, rndIdx);
    allResults(caseIdx).trueOrient = dipOrient;
    allResults(caseIdx).estLoc = estLoc;
    allResults(caseIdx).estOrient = estOrient;
    allResults(caseIdx).locError = locError;
    allResults(caseIdx).angleError = angleError;
    allResults(caseIdx).dipoleIndex = rndIdx;
end

%% ---- Part (خ): Compare Results for All Three Cases ----
fprintf('\n========================================\n');
fprintf('FINAL COMPARISON OF ALL THREE CASES\n');
fprintf('========================================\n');
fprintf('%-25s | %-15s | %-15s | %-15s\n', 'Case', 'Location Error (cm)', 'Angle Error (deg)', 'Dipole Index');
fprintf('%s\n', repmat('-', 75, 1));

for caseIdx = 1:length(allResults)
    fprintf('%-25s | %-15.4f | %-15.2f | %-15d\n', ...
        allResults(caseIdx).case, ...
        allResults(caseIdx).locError, ...
        allResults(caseIdx).angleError, ...
        allResults(caseIdx).dipoleIndex);
end

%% Create Comparison Figure for All Case
fig_comparison = figure('Name', 'Comparison of All Three Cases', 'Position', [50, 50, 1500, 500]);

for caseIdx = 1:length(allResults)
    subplot(1, 3, caseIdx);
    
    % Plot source grid
    scatter3(LocMat(1,:), LocMat(2,:), LocMat(3,:), 3, [0.2 0.6 1], 'filled');
    hold on;
    
    % Plot electrodes
    scatter3(ElecXYZ(:,1), ElecXYZ(:,2), ElecXYZ(:,3), 40, 'r', 'filled');
    
    % Plot true dipole
    trueLoc = allResults(caseIdx).trueLoc;
    trueOrient = allResults(caseIdx).trueOrient;
    quiver3(trueLoc(1), trueLoc(2), trueLoc(3), ...
            trueOrient(1), trueOrient(2), trueOrient(3), ...
            'g', 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    scatter3(trueLoc(1), trueLoc(2), trueLoc(3), 80, 'g', 'filled');
    
    % Plot estimated dipole
    estLoc = allResults(caseIdx).estLoc;
    estOrient = allResults(caseIdx).estOrient;
    quiver3(estLoc(1), estLoc(2), estLoc(3), ...
            estOrient(1), estOrient(2), estOrient(3), ...
            'm', 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    scatter3(estLoc(1), estLoc(2), estLoc(3), 80, 'm', 'filled');
    
    % Plot origin
    scatter3(0, 0, 0, 80, 'r', '*', 'LineWidth', 2);
    
    % Add error text
    textStr = sprintf('Loc Err: %.3f cm\nAng Err: %.1f°', ...
        allResults(caseIdx).locError, allResults(caseIdx).angleError);
    text(-11, -11, 11, textStr, 'BackgroundColor', 'w', 'FontSize', 9);
    
    xlabel('X (cm)'); ylabel('Y (cm)'); zlabel('Z (cm)');
    title(allResults(caseIdx).case);
    grid on; axis equal;
    view(45, 30);
    legend({'Source Grid', 'Electrodes', 'True', 'Est'}, 'Location', 'best', 'FontSize', 8);
end

sgtitle('Comparison of Source Localization for Three Dipole Configurations');

%% ---- Bar Plot of Errors ----
fig_errors = figure('Name', 'Error Comparison Bar Plot');
subplot(1,2,1);
bar([allResults.locError]);
set(gca, 'XTickLabel', {allResults.case});
ylabel('Location Error (cm)');
title('Localization Error Comparison');
grid on;

subplot(1,2,2);
bar([allResults.angleError]);
set(gca, 'XTickLabel', {allResults.case});
ylabel('Orientation Error (degrees)');
title('Orientation Error Comparison');
grid on;

%% ---- Save All Results as Images ----
% Create a folder to save images
saveFolder = 'EEG_Localization_Results';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Get all figure handles
allFigures = findall(0, 'Type', 'figure');

% Save each figure
for figIdx = 1:length(allFigures)
    figureHandle = allFigures(figIdx);
    
    % Get figure name
    figName = get(figureHandle, 'Name');
    if isempty(figName)
        figName = sprintf('Figure_%d', figIdx);
    end
    
    % Clean filename (remove invalid characters)
    figName = strrep(figName, ' ', '_');
    figName = strrep(figName, ':', '');
    figName = strrep(figName, '-', '_');
    
    % Save as PNG (high resolution)
    filename = fullfile(saveFolder, sprintf('%s.png', figName));
    saveas(figureHandle, filename);
    fprintf('Saved: %s\n', filename);
    
    % Also save as FIG (MATLAB format) for later editing
    filenameFig = fullfile(saveFolder, sprintf('%s.fig', figName));
    savefig(figureHandle, filenameFig);
end

% Save results data as MAT file
resultsFile = fullfile(saveFolder, 'localization_results.mat');
save(resultsFile, 'allResults', 'ModelParams', 'dipoleCases', 'caseNames');
fprintf('\nAll results saved to folder: %s\n', saveFolder);
fprintf('Results data saved to: %s\n', resultsFile);


%%  Bonus Section Optimization-based Dipole Localization 

% Prepare data for optimization (using the three dipole cases)
% We'll use the averaged potentials from Section 5 as measured data
optimResults = struct();
    % Get measured potentials for this case (use peakMean from earlier)

    % Note: We need to recompute peakMean for each case
    % For this example, we'll use the potentials we already have

for caseIdx = 1:length(allResults)
    % Get the true dipole information
    true_idx = allResults(caseIdx).dipoleIndex;
    true_loc = allResults(caseIdx).trueLoc;
    true_orient = allResults(caseIdx).trueOrient;
    
    % Generate measured potentials (simulated with noise optional)
    dip_cols_true = (true_idx-1)*3 + (1:3);
    G_true = GainMat(:, dip_cols_true);
    true_moment = 1e-6;  % Arbitrary magnitude
    measured_pot = G_true * (true_moment * true_orient);
    
    % Add small noise to make it realistic (2% noise)
    noise_level = 0.02 * max(abs(measured_pot(:)));
    measured_pot_noisy = measured_pot + noise_level * randn(size(measured_pot));
    
    %% Method 1: Simulated Annealing with 6 parameters (x,y,z,qx,qy,qz)
    fprintf('\n--- Method 1: Simulated Annealing (6 parameters) ---\n');
    
    % Bounds: [x, y, z, qx, qy, qz]
    r_max = ModelParams.R(1);
    q_max = 2e-6;  % Maximum dipole moment magnitude
    lb = [-r_max, -r_max, -r_max, -q_max, -q_max, -q_max];
    ub = [r_max, r_max, r_max, q_max, q_max, q_max];
    
    % Initial guess (random position near surface)
    init_params = [2*rand(1,3)-1; (2*rand(3,1)-1)*q_max]';
    init_params(1:3) = init_params(1:3) * r_max;
    
    % Simulated annealing options
    sa_options = optimoptions('simulannealbnd', ...
        'Display', 'iter', ...
        'MaxIterations', 500, ...
        'MaxStallIterations', 50, ...
        'TemperatureFcn', @temperatureexp, ...
        'AnnealingFcn', @annealingfast, ...
        'FunctionTolerance', 1e-6);
    
    % Run simulated annealing
    try
        [opt_params_sa, cost_sa] = simulannealbnd(@(p) dipole_cost_function_position_only(p, GainMat, LocMat, measured_pot_noisy, nElec), ...
            init_params, lb, ub, sa_options);
        
        % Extract optimized dipole
        opt_loc_sa = opt_params_sa(1:3)';
        opt_moment_sa = opt_params_sa(4:6)';
        opt_orient_sa = opt_moment_sa / norm(opt_moment_sa);
        
        % Find closest grid point
        distances_sa = sqrt(sum((LocMat - opt_loc_sa).^2, 1));
        [~, opt_idx_sa] = min(distances_sa);
        opt_loc_sa_grid = LocMat(:, opt_idx_sa);
        
        % Calculate errors
        loc_error_sa = norm(opt_loc_sa_grid - true_loc);
        angle_error_sa = acosd(abs(dot(opt_orient_sa, true_orient)));
        
        fprintf('SA Results:\n');
        fprintf('  Estimated location: [%.2f, %.2f, %.2f] cm\n', opt_loc_sa_grid(1), opt_loc_sa_grid(2), opt_loc_sa_grid(3));
        fprintf('  Location error: %.4f cm\n', loc_error_sa);
        fprintf('  Angle error: %.2f degrees\n', angle_error_sa);
        fprintf('  Final cost: %.6e\n', cost_sa);
        
    catch ME
        fprintf('Simulated Annealing failed: %s\n', ME.message);
        loc_error_sa = NaN;
        angle_error_sa = NaN;
        opt_loc_sa_grid = [NaN; NaN; NaN];
        opt_orient_sa = [NaN; NaN; NaN];
    end
    
end
    
    %% Store optimization results
    optimResults(caseIdx).case = allResults(caseIdx).case;
    optimResults(caseIdx).true_loc = true_loc;
    optimResults(caseIdx).true_orient = true_orient;
    
    optimResults(caseIdx).sa.loc = opt_loc_sa_grid;
    optimResults(caseIdx).sa.orient = opt_orient_sa;
    optimResults(caseIdx).sa.loc_error = loc_error_sa;
    optimResults(caseIdx).sa.angle_error = angle_error_sa;
    
    optimResults(caseIdx).ga.loc = opt_loc_ga;
    optimResults(caseIdx).ga.orient = opt_orient_ga;
    optimResults(caseIdx).ga.loc_error = loc_error_ga;
    optimResults(caseIdx).ga.angle_error = angle_error_ga;
    
    optimResults(caseIdx).ps.loc = opt_loc_ps_grid;
    optimResults(caseIdx).ps.orient = opt_orient_ps;
    optimResults(caseIdx).ps.loc_error = loc_error_ps;
    optimResults(caseIdx).ps.angle_error = angle_error_ps;
    
    %% Visualize optimization results for this case
    fig_optim = figure('Name', sprintf('Optimization Results - %s', allResults(caseIdx).case));
    scatter3(LocMat(1,:), LocMat(2,:), LocMat(3,:), 5, [0.2 0.6 1], 'filled');
    hold on;
    
    % Plot electrodes
    for k = 1:nElec
        scatter3(ElecXYZ(k,1), ElecXYZ(k,2), ElecXYZ(k,3), 40, 'r', 'filled');
    end
    
    % Plot true dipole
    quiver3(true_loc(1), true_loc(2), true_loc(3), ...
            true_orient(1), true_orient(2), true_orient(3), ...
            'g', 'LineWidth', 2.5, 'MaxHeadSize', 1.5);
    
    % Plot SA result
    if ~isnan(loc_error_sa)
        quiver3(opt_loc_sa_grid(1), opt_loc_sa_grid(2), opt_loc_sa_grid(3), ...
                opt_orient_sa(1), opt_orient_sa(2), opt_orient_sa(3), ...
                'm', 'LineWidth', 2, 'MaxHeadSize', 1.5, 'LineStyle', '-');
    end
    
    % Plot GA result
    if ~isnan(loc_error_ga)
        quiver3(opt_loc_ga(1), opt_loc_ga(2), opt_loc_ga(3), ...
                opt_orient_ga(1), opt_orient_ga(2), opt_orient_ga(3), ...
                'c', 'LineWidth', 2, 'MaxHeadSize', 1.5, 'LineStyle', '--');
    end
    
    % Plot Pattern Search result
    if ~isnan(loc_error_ps)
        quiver3(opt_loc_ps_grid(1), opt_loc_ps_grid(2), opt_loc_ps_grid(3), ...
                opt_orient_ps(1), opt_orient_ps(2), opt_orient_ps(3), ...
                'y', 'LineWidth', 2, 'MaxHeadSize', 1.5, 'LineStyle', ':');
    end
    
    scatter3(0, 0, 0, 80, 'r', '*', 'LineWidth', 2);
    
    xlabel('X (cm)'); ylabel('Y (cm)'); zlabel('Z (cm)');
    title(sprintf('Optimization Results - %s', allResults(caseIdx).case));
    legend({'Source Grid', 'Electrodes', 'True Dipole', 'Simulated Annealing', 'Genetic Algorithm', 'Pattern Search', 'Origin'}, ...
           'Location', 'best');
    grid on; axis equal;
    view(45, 30);
    

% Save optimization results
optimResultsFile = fullfile(saveFolder, 'optimization_results.mat');
save(optimResultsFile, 'optimResults', 'comparisonTable');
fprintf('\nOptimization results saved to: %s\n', optimResultsFile);

% Save comparison figure
saveas(fig_comparison_optim, fullfile(saveFolder, 'MNE_vs_Optimization_Comparison.png'));
fprintf('Comparison figure saved\n');

fprintf('\n========================================\n');
fprintf('BONUS SECTION COMPLETED\n');

%% Define the cost function for optimization
% Cost function: ||m - G*q||^2
% m: measured potentials (nElec x 1)
% G: gain matrix for a specific dipole position (nElec x 3)
% q: dipole moment vector (3 x 1)

function cost = dipole_cost_function_position_only(params, GainMat, LocMat, measured_pot, nElec)
    % params: [dipole_index, qx, qy, qz]
    % OR: params: [x, y, z, qx, qy, qz]
    
    if length(params) == 4
        % Case: dipole index + moment
        dipole_idx = round(params(1));
        dipole_idx = max(1, min(dipole_idx, size(LocMat, 2)));
        q = params(2:4)';
    else
        % Case: x,y,z + moment
        x = params(1); y = params(2); z = params(3);
        dipole_loc = [x; y; z];
        distances = sqrt(sum((LocMat - dipole_loc).^2, 1));
        [~, dipole_idx] = min(distances);
        q = params(4:6)';
    end
    
    % Get lead field for this dipole
    dip_cols = (dipole_idx-1)*3 + (1:3);
    G_dip = GainMat(:, dip_cols);
    
    % Compute forward potentials
    computed_pot = G_dip * q;
    
    % Cost = squared error
    cost = norm(measured_pot - computed_pot)^2;
end

% Alternative: Cost function with normalized orientation (5 parameters)
function cost = dipole_cost_function_orientation(params, GainMat, LocMat, measured_pot, nElec)
    % params: [dipole_idx, theta, phi, magnitude]
    % theta: polar angle (0 to pi), phi: azimuthal angle (0 to 2pi)
    
    dipole_idx = round(params(1));
    dipole_idx = max(1, min(dipole_idx, size(LocMat, 2)));
    theta = params(2);
    phi = params(3);
    magnitude = abs(params(4));
    
    % Orientation vector
    orient = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
    orient = orient / norm(orient);
    
    % Dipole moment
    q = magnitude * orient;
    
    % Get lead field
    dip_cols = (dipole_idx-1)*3 + (1:3);
    G_dip = GainMat(:, dip_cols);
    
    % Compute forward potentials
    computed_pot = G_dip * q;
    
    % Cost = squared error
    cost = norm(measured_pot - computed_pot)^2;
end

