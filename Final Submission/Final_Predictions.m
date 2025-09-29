%% Time Setup
startTime = datetime(2025, 9, 26, 2, 1, 6, 744000);
stopTime = startTime + hours(48);
sampleTime = 100; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

%% Earth Parameters
earthRadius = 6.378137e6; % meters
mu_earth = 3.986004418e14;

N = 5;
numSamples = 10000;
lat = 10;
lon = -40;

sc_mean = satelliteScenario(startTime, stopTime, sampleTime);

% --------------------------
% Processing Object 1 state
% --------------------------

fprintf('\n=== Processing Object1 ===\n');

Object1 = 'CZ4C DEB.xml';

sat = satellite(sc_mean, Object1);
elements1 = orbitalElements(sat);

% Use elements1.Epoch for Object1 timing
startTime_Object1 = elements1.Epoch;
stopTime_Object1 = startTime_Object1 + days(2);
timeVector1 = startTime_Object1:seconds(sampleTime):stopTime_Object1;
numTimeSteps_Object1 = length(timeVector1);  % Use the correct length

[nominal, P_Object1] = classical_to_equinoctial(elements1);

mu_Object1 = [nominal.a;
nominal.h;
nominal.k;
nominal.p;
nominal.q;
nominal.l];

[Object1_weighted_covariance, Object1_weighted_mean, Object1_gaussian_means, Object1_gaussian_covariances, ~, ~, Object1_gaussian_weights, ~, ~] = Generate_GSUKF(mu_Object1, P_Object1, N, numSamples, startTime_Object1, stopTime_Object1, sampleTime, Object1);

Object1_positions = Object1_weighted_mean(1:3, :);
Object1_position_covs = Object1_weighted_covariance(1:3, 1:3, :);
Object1_positions_km = Object1_positions / 1000;
time_hours_Object1 = hours(timeVector1 - startTime_Object1);

% ----------------------------
% Plotting Object 1 positions
% ----------------------------

figure('Position', [200, 200, 1400, 1000]);

position_magnitude_Object1 = sqrt(sum(Object1_positions_km.^2, 1));
fprintf('Object1 Position magnitude size: %dx%d\n', size(position_magnitude_Object1));
fprintf('Time hours size: %dx%d\n', size(time_hours_Object1));
fprintf('Number of time steps for Object1: %d\n', numTimeSteps_Object1);

plot(timeVector1, position_magnitude_Object1, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Object1 Position Magnitude');
hold on;

% Calculate uncertainty bounds for Object1 - use correct array size
std_dev_magnitude_Object1 = zeros(1, numTimeSteps_Object1);  % Fixed size
for i = 1:numTimeSteps_Object1
    combined_variance = sum(diag(Object1_position_covs(:, :, i))) / (1000^2);
    std_dev_magnitude_Object1(i) = sqrt(combined_variance);
end

% Add uncertainty bounds for Object1
upper_bound_Object1 = position_magnitude_Object1 + std_dev_magnitude_Object1;
lower_bound_Object1 = position_magnitude_Object1 - std_dev_magnitude_Object1;
fill([timeVector1, fliplr(timeVector1)], [upper_bound_Object1, fliplr(lower_bound_Object1)], ...
    'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('Time (hours)');
ylabel('Position Magnitude (km)');
title('Object1 Position Magnitude vs Time (±1σ)');
legend('Location', 'best');
grid on;

% ----------------------------------------------------
% Plotting Object 1 Gaussians and Monte Carlo samples
% ----------------------------------------------------

numSamples_MC = 300;

[next_Particle_ECI, ~] = Generate_Monte_Carlo(mu_Object1, P_Object1, startTime_Object1, stopTime_Object1, sampleTime, numSamples_MC, Object1);

% Convert positions to km
Positions_MC_km = next_Particle_ECI(1:3, :, :) / 1000;
GSF_positions_km = Object1_weighted_mean(1:3, :) / 1000;
GSF_position_covs = Object1_weighted_covariance(1:3, 1:3, :);

% Calculate Monte Carlo statistics
MC_mean_positions = squeeze(mean(Positions_MC_km, 3));
MC_std_positions = squeeze(std(Positions_MC_km, 0, 3));

% Calculate position magnitudes
MC_magnitudes = squeeze(sqrt(sum(Positions_MC_km.^2, 1)));
MC_mean_magnitude = mean(MC_magnitudes, 2);
MC_std_magnitude = std(MC_magnitudes, 0, 2);

% GSF magnitude calculations
GSF_magnitude = sqrt(sum(GSF_positions_km.^2, 1));
GSF_magnitude_std = zeros(size(GSF_magnitude));
for t = 1:size(GSF_position_covs, 3)
    GSF_magnitude_std(t) = sqrt(trace(GSF_position_covs(:,:,t))) / (1000);
end

% Ensure consistent dimensions
min_length = min([length(MC_mean_magnitude), length(GSF_magnitude)]);
MC_mean_magnitude = MC_mean_magnitude(1:min_length);
MC_std_magnitude = MC_std_magnitude(1:min_length);
GSF_magnitude = GSF_magnitude(1:min_length);
GSF_magnitude_std = GSF_magnitude_std(1:min_length);

% Time vector for plotting (ensure row vector)
time_hours = hours(timeVector1 - startTime_Object1);
time_hours = time_hours(:)'; % Force row vector
time_hours = time_hours(1:min_length); % Match other arrays

%% ======== VISUALIZATION ==========
figure('Position', [200, 200, 1400, 1000]);
analysis_time_hours = 4;
analysis_index = round(analysis_time_hours * 3600 / sampleTime) + 1;

if analysis_index <= numTimeSteps
    
    % Monte Carlo samples at this time
    MC_samples_analysis = squeeze(Positions_MC_km(:, analysis_index, :))'; 
    
    % Plot Monte Carlo samples
    scatter3(MC_samples_analysis(:,1), MC_samples_analysis(:,2), MC_samples_analysis(:,3), ...
        10, 'b', 'filled', 'MarkerFaceAlpha', 0.2, 'DisplayName', 'MC Samples');
    hold on;
    
    % Plot each individual Gaussian component
    component_colors = lines(N);
    for comp = 1:size(Object1_gaussian_means, 3)
        comp_mean = Object1_gaussian_means(1:3, analysis_index, comp) / 1000;
        comp_cov = Object1_gaussian_covariances(1:3, 1:3, analysis_index, comp) / (1000^2);
        comp_weight = Object1_gaussian_weights(comp);
        
        % Plot component mean with size proportional to weight
        scatter3(comp_mean(1), comp_mean(2), comp_mean(3), ...
            100 + 200*comp_weight, component_colors(comp,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('G%d ', comp));
        
        % Plot 1-sigma ellipsoid with transparency proportional to weight
        try
            plot_uncertainty_ellipsoid(comp_mean, comp_cov, 1, component_colors(comp,:), ...
                0.1 + 0.4*comp_weight);
        catch ME
            fprintf('Warning: Could not plot ellipsoid for component %d: %s\n', comp, ME.message);
        end
    end
    
    title(sprintf('GSF Components at t = %d hours', analysis_time_hours));
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
    legend('Location', 'best', 'FontSize', 8);
    grid on; axis equal;
end

% -------------------------
% Processing Object2 state
% -------------------------

fprintf('\n=== Processing Object2 ===\n');

Object2 = 'COSMOS 1570.xml';

sat = satellite(sc_mean, Object2);
elements2 = orbitalElements(sat);

% Use elements2.Epoch for Object2 timing
startTime_Object2 = elements2.Epoch;
stopTime_Object2 = startTime_Object2 + days(2);
timeVector2 = startTime_Object2:seconds(sampleTime):stopTime_Object2;
numTimeSteps_Object2 = length(timeVector2);  % Use the correct length

[nominal, P_Object2] = classical_to_equinoctial(elements2);

mu_Object2 = [nominal.a;
    nominal.h;
    nominal.k;
    nominal.p;
    nominal.q;
    nominal.l];

[Object2_weighted_covariance, Object2_weighted_mean, Object2_gaussian_means, Object2_gaussian_covariances, ~, ~, Object2_gaussian_weights, ~, ~] = Generate_GSUKF(mu_Object2, P_Object2, N, numSamples, startTime_Object2, stopTime_Object2, sampleTime, Object2);

Object2_positions = Object2_weighted_mean(1:3, :);
Object2_position_covs = Object2_weighted_covariance(1:3, 1:3, :);
Object2_positions_km = Object2_positions / 1000;
time_hours_Object2 = hours(timeVector2 - startTime_Object2);

% ---------------------------
% Plotting Object2 positions
% ---------------------------

figure('Position', [200, 200, 1400, 1000]);

position_magnitude_Object2 = sqrt(sum(Object2_positions_km.^2, 1));
fprintf('Object2 Position magnitude size: %dx%d\n', size(position_magnitude_Object2));
fprintf('Number of time steps for Object2: %d\n', numTimeSteps_Object2);

plot(timeVector2, position_magnitude_Object2, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Object2 Position Magnitude');
hold on;

% Calculate uncertainty bounds for Object2 - use correct array size
std_dev_magnitude_Object2 = zeros(1, numTimeSteps_Object2);  % Fixed size
for i = 1:numTimeSteps_Object2
    combined_variance = sum(diag(Object2_position_covs(:, :, i))) / (1000^2);
    std_dev_magnitude_Object2(i) = sqrt(combined_variance);
end

% Add uncertainty bounds for Object2
upper_bound_Object2 = position_magnitude_Object2 + std_dev_magnitude_Object2;
lower_bound_Object2 = position_magnitude_Object2 - std_dev_magnitude_Object2;
fill([timeVector2, fliplr(timeVector2)], [upper_bound_Object2, fliplr(lower_bound_Object2)], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

xlabel('Time (hours)');
ylabel('Position Magnitude (km)');
title('Object2 Position Magnitude vs Time (±1σ)');
legend('Location', 'best');
grid on;

% ---------------------------------
% Finding overlaps between objects
% ---------------------------------

fprintf('\n=== Calculating Relative Distance ===\n');

% Find the overlapping time period
common_start = max(timeVector1(1), timeVector2(1));
common_end = min(timeVector1(end), timeVector2(end));

fprintf('Object1 time range: %s to %s\n', datestr(timeVector1(1)), datestr(timeVector1(end)));
fprintf('Object2 time range: %s to %s\n', datestr(timeVector2(1)), datestr(timeVector2(end)));
fprintf('Common time range: %s to %s\n', datestr(common_start), datestr(common_end));

% Check if there's actually an overlap
if common_start >= common_end
    error('No time overlap between the two satellites!');
end

% Create a common time vector for the overlapping period
common_timeVector = common_start:seconds(sampleTime):common_end;
numCommonSteps = length(common_timeVector);

fprintf('Number of common time steps: %d\n', numCommonSteps);

% For Object1
Object1_common = zeros(6, numCommonSteps);
Object1_common_covs = zeros(6, 6, numCommonSteps);

for i = 1:6
    Object1_common(i, :) = interp1(timeVector1, Object1_weighted_mean(i, :), common_timeVector, 'linear', 'extrap');
end

for i = 1:6
    for j = 1:6
        cov_series = squeeze(Object1_weighted_covariance(i, j, :));
        Object1_common_covs(i, j, :) = interp1(timeVector1, cov_series, common_timeVector, 'linear', 'extrap');
    end
end

% For Object2
Object2_common = zeros(6, numCommonSteps);
Object2_common_covs = zeros(6, 6, numCommonSteps);

for i = 1:6
    Object2_common(i, :) = interp1(timeVector2, Object2_weighted_mean(i, :), common_timeVector, 'linear', 'extrap');
end

for i = 1:3
    for j = 1:3
        cov_series = squeeze(Object2_weighted_covariance(i, j, :));
        Object2_common_covs(i, j, :) = interp1(timeVector2, cov_series, common_timeVector, 'linear', 'extrap');
    end
end


%% Interpolate Object1 Gaussian Components
Object1_gaussian_means_common = zeros(6, numCommonSteps, N);
Object1_gaussian_covariances_common = zeros(6, 6, numCommonSteps, N);

% Interpolate means for each Gaussian component
for g = 1:N
    for i = 1:6
        mean_series = squeeze(Object1_gaussian_means(i, :, g));
        Object1_gaussian_means_common(i, :, g) = interp1(timeVector1, mean_series, common_timeVector, 'linear', 'extrap');
    end
end

% Interpolate covariances for each Gaussian component
for g = 1:N
    for i = 1:6
        for j = 1:6
            cov_series = squeeze(Object1_gaussian_covariances(i, j, :, g));
            Object1_gaussian_covariances_common(i, j, :, g) = interp1(timeVector1, cov_series, common_timeVector, 'linear', 'extrap');
        end
    end
end


%% Interpolate Object2 Gaussian Components
Object2_gaussian_means_common = zeros(6, numCommonSteps, N);
Object2_gaussian_covariances_common = zeros(6, 6, numCommonSteps, N);

% Interpolate means for each Gaussian component
for g = 1:N
    for i = 1:6
        mean_series = squeeze(Object2_gaussian_means(i, :, g));
        Object2_gaussian_means_common(i, :, g) = interp1(timeVector2, mean_series, common_timeVector, 'linear', 'extrap');
    end
end

% Interpolate covariances for each Gaussian component
for g = 1:N
    for i = 1:6
        for j = 1:6
            cov_series = squeeze(Object2_gaussian_covariances(i, j, :, g));
            Object2_gaussian_covariances_common(i, j, :, g) = interp1(timeVector2, cov_series, common_timeVector, 'linear', 'extrap');
        end
    end
end

%% Calculate Relative Distance and Uncertainty
relative_positions = Object1_common(1:3, :) - Object2_common(1:3, :);
relative_distances = sqrt(sum(relative_positions.^2, 1)) / 1000; % Convert to km

%% Plot Results
figure('Position', [300, 300, 1400, 1000]);

% Plot relative distance
plot(common_timeVector, relative_distances, 'r-', 'LineWidth', 2, 'DisplayName', 'Relative Distance');
hold on;

% Add some statistics
fprintf('\nRelative Distance Statistics:\n');
fprintf('Minimum distance: %.2f km\n', min(relative_distances));
fprintf('Maximum distance: %.2f km\n', max(relative_distances));
fprintf('Mean distance: %.2f km\n', mean(relative_distances));
fprintf('Standard deviation of distance: %.2f km\n', std(relative_distances));

% Find closest approach
[min_distance, min_idx] = min(relative_distances);
closest_time = common_timeVector(min_idx);
fprintf('Closest approach: %.2f km at %s\n', min_distance, datestr(closest_time));

R = 6; % Combined hard body radius of the two objects.

% Call the modified propagate function
[collisionprobs_total, lognoprob, noprob] = propagate_gaussian_mixture(... 
    Object1_gaussian_means_common, Object2_gaussian_means_common , ...
    Object1_gaussian_covariances_common , Object2_gaussian_covariances_common , ...
    Object1_gaussian_weights, Object2_gaussian_weights, ...
    numCommonSteps, N, R);

fprintf('Calculating total collision probability across all debris...\n');
for j = 1:numCommonSteps
    lognoprob(j) = 0; 
    for i = 1:numDebris
        lognoprob(j) = lognoprob(j) + log1p(-collisionprobs_total(1, j, i));
    end
    noprob(j) = exp(lognoprob(j));
end
total_collision_prob_per_timestep = 1 - noprob; 
max_total_prob = max(total_collision_prob_per_timestep);

fprintf('\n=== TOTAL COLLISION PROBABILITY ANALYSIS ===\n');
fprintf('Maximum total collision probability: %.6e\n', max_total_prob);

[~, max_idx] = max(total_collision_prob_per_timestep);
fprintf('Maximum occurs at time step: %d\n', max_idx);
fprintf('Time: %s\n', char(common_timeVector(max_idx)));

high_risk_total = sum(total_collision_prob_per_timestep > 1e-4);
if high_risk_total > 0
    fprintf('HIGH RISK: %d time steps with total Pc > 10^-4\n', high_risk_total);
end

% Plotting code remains the same...
figure('Position', [100, 100, 1200, 600]);
plot(common_timeVector, total_collision_prob_per_timestep, 'b-', 'LineWidth', 2);
hold on;

yline(1e-4, 'm--', 'LineWidth', 1.5, 'Label', 'Pc = 10^{-4}');

[max_total_prob, max_idx] = max(total_collision_prob_per_timestep);
plot(common_timeVector(max_idx), max_total_prob, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

text(common_timeVector(max_idx), max_total_prob, ...
     sprintf('Max Pc = %.2e', max_total_prob), ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10);

xlabel('Time');
ylabel('Total Collision Probability');
title('Total Collision Probability over Time - Linear Scale');
grid on;
hold off;


