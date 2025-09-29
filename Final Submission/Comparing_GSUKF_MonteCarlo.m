%% Time Setup
startTime = datetime(2025, 9, 26, 2, 1, 6, 744000);
stopTime = startTime + hours(24);
sampleTime = 800; % seconds
timeVector = startTime:seconds(sampleTime):stopTime;
numTimeSteps = length(timeVector);

%% Earth Parameters
earthRadius = 6.378137e6; % meters
mu_earth = 3.986004418e14;

numSamples = 10000;
num_Monte_Carlo_samples = 1000;
lat = 10;
lon = -50;

sc_mean = satelliteScenario(startTime, stopTime, sampleTime);

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

% ----------------------------------------------
% Generating the ground truth Monte Carlo state
% ----------------------------------------------

tic_total = tic;

[next_Particle_ECI, Weights] = Generate_Monte_Carlo(mu_Object1, P_Object1, startTime_Object1, stopTime_Object1, sampleTime, num_Monte_Carlo_samples, Object1);

% Initialize arrays for statistical calculations
Particle_weighted_mean = zeros(6, numTimeSteps_Object1);
Particle_weighted_covariance = zeros(6, 6, numTimeSteps_Object1);
Particle_weighted_skewness = zeros(6, numTimeSteps_Object1);
Particle_weighted_kurtosis = zeros(6, numTimeSteps_Object1);

% Statistical calculations loop - corrected version
for t = 1:numTimeSteps_Object1
    % Extract particles at time step t: 6 x numParticles
    particles_t = squeeze(next_Particle_ECI(:, t, :)); % 6 x numParticles
    
    % Weights normalized and in correct shape
    w = ones(1, num_Monte_Carlo_samples) / num_Monte_Carlo_samples;

    % Weighted mean: 6 x 1 vector
    Particle_weighted_mean(:, t) = particles_t * w';
    
    % Center the particles
    particles_centered = particles_t - Particle_weighted_mean(:, t);
    
    % For weighted covariance with bias correction
    % Using the standard weighted covariance formula
    Particle_weighted_covariance(:, :, t) = (particles_centered .* w) * particles_centered' / (1 - sum(w.^2));
    
    % Weighted skewness and kurtosis for each component
    for j = 1:6
        x = particles_centered(j, :); % 1 x numParticles - already centered!
        
        % Calculate weighted central moments
        W1 = sum(w);           % Sum of weights (should be 1)
        W2 = sum(w.^2);        % Sum of squared weights
        W3 = sum(w.^3);        % Sum of cubed weights
        
        % Weighted central moments
        mu2 = sum(w .* (x.^2)) / W1;  % 2nd central moment
        mu3 = sum(w .* (x.^3)) / W1;  % 3rd central moment  
        mu4 = sum(w .* (x.^4)) / W1;  % 4th central moment
        
        % Check for near-zero variance to avoid division by zero
        if mu2 < eps
            Particle_weighted_skewness(j, t) = 0;
            Particle_weighted_kurtosis(j, t) = 3; % Normal distribution kurtosis
        else
            % Weighted skewness
            Particle_weighted_skewness(j, t) = mu3 / (mu2^(3/2));
            
            % Weighted kurtosis (raw kurtosis)
            Particle_weighted_kurtosis(j, t) = mu4 / (mu2^2);
        end
    end
end

% Extract mean positions and covariances
Particle_ISS_positions = Particle_weighted_mean(1:3, :);
Particle_ISS_position_covs = Particle_weighted_covariance(1:3, 1:3, :);
Particle_ISS_positions_skewness = Particle_weighted_skewness(1:3, :);
Particle_ISS_positions_kurtosis = Particle_weighted_kurtosis(1:3, :);
variances = squeeze([Particle_ISS_position_covs(1,1,end); ...
                     Particle_ISS_position_covs(2,2,end); ...
                     Particle_ISS_position_covs(3,3,end)]);

% Convert to uncertainties (standard deviations)
uncertainties = sqrt(variances);

% Uncertainty in the MEAN (standard error of the mean)
N = size(Particle_weighted_mean, 2);  % number of MC samples
mean_cov = Particle_ISS_position_covs(:,:,end) / N;   % covariance of the mean
mean_std = sqrt(diag(mean_cov));      % SE per coordinate

disp('Mean position at final time:');
disp(Particle_ISS_positions(:, end));

disp('1-sigma ensemble spread (m):');
disp(uncertainties);

disp('1-sigma uncertainty in the mean (standard error, m):');
disp(mean_std);

%% Display total computation time
total_time = toc(tic_total);
fprintf('\nTotal computation time: %.2f seconds (%.2f minutes)\n', total_time, total_time/60);


%% Define different numbers of Gaussian components to test
N_values = [1, 20, 40, 60, 80];

fprintf('\n=== Processing Multiple N Values ===\n');

% Store results for each N value
GSF_results = cell(length(N_values), 1);
bhattacharyya_distances = cell(length(N_values), 1);
euclidean_distances = cell(length(N_values), 1);
% NEW: Store skewness and kurtosis results
GSF_skewness = cell(length(N_values), 1);
GSF_kurtosis = cell(length(N_values), 1);
computation_times = zeros(length(N_values), 1);

% Generate particle filter reference (or use the highest N as reference)
ref_positions = Particle_ISS_positions(1:3, :);
ref_position_covs = Particle_ISS_position_covs(1:3, 1:3, :);
ref_position_skewness = Particle_ISS_positions_skewness(1:3, :);
ref_position_kurtosis = Particle_ISS_positions_kurtosis(1:3, :);

% Process each N value
for idx = 1:length(N_values)
    N = N_values(idx);
    fprintf('Processing N = %d Gaussians...\n', N);
    
    % Start timing
    tic;
    
    [weighted_covariance, weighted_mean, gaussian_ECI_means, gaussian_ECI_covs, ~, ~, Weights, Wm_all, Wc_all] = Generate_GSUKF(mu_Object1, P_Object1, N, numSamples, startTime_Object1, stopTime_Object1, sampleTime, Object1);
   
    % End timing
    computation_times(idx) = toc;
    
    fprintf('  Computation time for N=%d: %.3f seconds\n', N, computation_times(idx));
    
    timeVector = datetime(timeVector, 'TimeZone', '');

    % Store results
    GSF_results{idx}.weighted_covariance = weighted_covariance;
    GSF_results{idx}.weighted_mean = weighted_mean;
    GSF_results{idx}.N = N;
    
    % Extract position data
    ISS_positions = weighted_mean(1:3, :);
    ISS_position_covs = weighted_covariance(1:3, 1:3, :);
    
    % NEW: Calculate skewness and kurtosis for GS-UKF
    [gsf_skewness, gsf_kurtosis] = calculate_mixture_moments(gaussian_ECI_means(1:3,:,:), gaussian_ECI_covs(1:3,1:3,:,:), Weights, numTimeSteps_Object1);
    GSF_skewness{idx} = gsf_skewness;
    GSF_kurtosis{idx} = gsf_kurtosis;
    
    % Calculate Bhattacharyya distance against reference
    d_bhatt = zeros(1, numTimeSteps_Object1);
    d_eucl = zeros(1, numTimeSteps_Object1);
    
    for k = 1:numTimeSteps_Object1
        % Bhattacharyya distance calculation
        mu1 = ISS_positions(:, k);
        Sigma1 = ISS_position_covs(:, :, k);
        mu2 = ref_positions(:, k);
        Sigma2 = ref_position_covs(:, :, k);
        
        d_bhatt(k) = bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2);
        
        % Euclidean distance calculation
        d_eucl(k) = norm(ISS_positions(:, k) - ref_positions(:, k));
    end
    
    bhattacharyya_distances{idx} = d_bhatt;
    euclidean_distances{idx} = d_eucl;
    
    fprintf('  Mean Bhattacharyya distance: %.6f\n', mean(d_bhatt));
    fprintf('  Max Bhattacharyya distance: %.6f\n', max(d_bhatt));
    fprintf('  Mean Euclidean distance: %.3f m\n', mean(d_eucl));
    fprintf('  Max Euclidean distance: %.3f m\n', max(d_eucl));
    fprintf('  RMS Euclidean error: %.3f m\n', sqrt(mean(d_eucl.^2)));
    
    % NEW: Print skewness and kurtosis statistics
    fprintf('  Mean position skewness (X,Y,Z): [%.3f, %.3f, %.3f]\n', ...
        mean(gsf_skewness(1,:)), mean(gsf_skewness(2,:)), mean(gsf_skewness(3,:)));
    fprintf('  Mean position kurtosis (X,Y,Z): [%.3f, %.3f, %.3f]\n', ...
        mean(gsf_kurtosis(1,:)), mean(gsf_kurtosis(2,:)), mean(gsf_kurtosis(3,:)));
end

%% Display computation time statistics
fprintf('\n=== Computation Time Analysis ===\n');
fprintf('%-8s %-15s\n', 'N', 'Time (seconds)');
fprintf('%-8s %-15s\n', '---', '-------------');
for idx = 1:length(N_values)
    fprintf('%-8d %-15.3f\n', N_values(idx), computation_times(idx));
end

%% Plotting
%timeVector = datetime(timeVector, 'TimeZone', '');
time_hours = hours(timeVector1 - startTime_Object1);

% Convert reference positions to km for plotting
ref_positions_km = ref_positions / 1000;
earthRadius_km = earthRadius / 1000;

figure('Position', [100, 100, 1400, 1000]);
% Plot 1: Mean Bhattacharyya Distance vs N
subplot(2, 2, 1);
mean_distances = cellfun(@mean, bhattacharyya_distances);
semilogx(N_values, mean_distances, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Mean Bhattacharyya Distance');
title('GS-UKF Convergence: Mean Distance vs N');
grid on;

% Plot 2: Max Bhattacharyya Distance vs N  
subplot(2, 2, 2);
max_distances = cellfun(@max, bhattacharyya_distances);
semilogx(N_values, max_distances, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Maximum Bhattacharyya Distance');
title('GS-UKF Convergence: Maximum Distance vs N');
grid on;

% Plot 3: Standard deviation of distances vs N
subplot(2, 2, 3);
std_distances = cellfun(@std, bhattacharyya_distances);
semilogx(N_values, std_distances, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Std Dev of Bhattacharyya Distance');
title('GS-UKF Variability: Std Dev vs N');
grid on;

% Plot 4: Distance evolution over time for all N values
subplot(2, 2, 4);
hold on;
for idx = 1:length(N_values)
    d = bhattacharyya_distances{idx};
    plot(time_hours, d, ...
         'LineWidth', 1.5, 'DisplayName', sprintf('N=%d', N_values(idx)));
end
xlabel('Time (hours)');
ylabel('Bhattacharyya Distance');
title('GS-UKF Distance Evolution: All N Values');
legend('Location', 'best');
grid on;


figure('Position', [150, 150, 1400, 1000]);
% Plot 1: Mean Euclidean Distance vs N
subplot(2, 2, 1);
mean_distances = cellfun(@mean, euclidean_distances);
semilogx(N_values, mean_distances, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Mean Euclidean Distance (m)');
title('GS-UKF Convergence: Mean Euclidean Distance vs N');
grid on;

% Plot 2: Max Euclidean Distance vs N
subplot(2, 2, 2);
max_distances = cellfun(@max, euclidean_distances);
semilogx(N_values, max_distances, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Max Euclidean Distance (m)');
title('GS-UKF Convergence: Max Euclidean Distance vs N');
grid on;

% Plot 3: Std Dev of Euclidean Distance vs N
subplot(2, 2, 3);
std_distances = cellfun(@std, euclidean_distances);
semilogx(N_values, std_distances, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Gaussian Components (N)');
ylabel('Std Dev of Euclidean Distance (m)');
title('GS-UKF Variability: Std Dev vs N');
grid on;

% Plot 4: Euclidean Distance Evolution Over Time
subplot(2, 2, 4);
hold on;
for idx = 1:length(N_values)
    d = euclidean_distances{idx};
    plot(time_hours, d, 'LineWidth', 1.5, 'DisplayName', sprintf('N=%d', N_values(idx)));
end
xlabel('Time (hours)');
ylabel('Euclidean Distance (m)');
title('GS-UKF Distance Evolution: All N Values');
legend('Location', 'best');
grid on;


figure('Position', [50, 50, 1400, 900]);
% Computation time with different N
plot(N_values, computation_times, 'bo-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
xlabel('Number of Gaussian Components (N)');
ylabel('Computation Time (seconds)');
title('Computation Time vs Number of Gaussians');
grid on;


figure('Position', [250, 250, 1600, 1200]);
% Plot 1: Skewness
subplot(2, 1, 1);
hold on;
for idx = 1:length(N_values)
    skew = GSF_skewness{idx}(3,:);
    plot(time_hours, skew, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('N=%d', N_values(idx)));
end
plot(time_hours, ref_position_skewness(3,:), 'k--', 'Color', [0 0 0 0.2], 'LineWidth', 0.2, ...
     'DisplayName', 'Reference');
xlabel('Time (hours)');
ylabel('Skewness (X-component)');
title('Skewness Evolution: X-component');
legend('Location', 'best');
grid on;

subplot(2, 1, 2);
% Plot 2: Kurtosis
hold on;
for idx = 1:length(N_values)
    kurt = GSF_kurtosis{idx}(3,:);
    plot(time_hours, kurt, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('N=%d', N_values(idx)));
end
plot(time_hours, ref_position_kurtosis(3,:), 'k--', 'Color', [0 0 0 0.2], 'LineWidth', 1, ...
     'DisplayName', 'Reference');
xlabel('Time (hours)');
ylabel('Kurtosis (X-component)');
title('Kurtosis Evolution: X-component');
legend('Location', 'best');
grid on;


figure('Position', [350, 350, 1400, 1000]);
% Covariance volume with different N
hold on;
ref_det = zeros(numTimeSteps_Object1, 1);
for t = 1:numTimeSteps_Object1
    ref_det(t) = det(ref_position_covs(:,:,t));
end
semilogy(time_hours, ref_det.^(1/6) / 1000, 'k--', 'LineWidth', 2, 'DisplayName', 'Reference');
for idx = 1:length(N_values)
    gsf_cov = GSF_results{idx}.weighted_covariance(1:3, 1:3, :);
    gsf_det = zeros(numTimeSteps_Object1, 1);
    for t = 1:numTimeSteps_Object1
        gsf_det(t) = det(gsf_cov(:,:,t));
    end
    semilogy(time_hours, gsf_det.^(1/6) / 1000, 'LineWidth', 1.5, ...
             'DisplayName', sprintf('N=%d', N_values(idx)));
end
xlabel('Time (hours)');
ylabel('Uncertainty Volume^{1/6} (km)');
title('GS-UKF Covariance Determinant Evolution');
legend('Location', 'best');
grid on;

for idx = 1:length(N_values)
    d_eucl = euclidean_distances{idx};
    gsf_skew = GSF_skewness{idx};
    gsf_kurt = GSF_kurtosis{idx};
    
    eucl_rms = sqrt(mean(d_eucl.^2));
    skew_rms = sqrt(mean((gsf_skew - ref_position_skewness).^2, 'all'));
    kurt_rms = sqrt(mean((gsf_kurt - ref_position_kurtosis).^2, 'all'));
    
    fprintf('%-8d %-12.3f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-15.3f\n', ...
        N_values(idx), eucl_rms, skew_rms, kurt_rms, ...
        mean(gsf_skew(1,:)), mean(gsf_skew(2,:)), mean(gsf_skew(3,:)), ...
        computation_times(idx));
end

fprintf('\nReference Values:\n');
fprintf('Skewness (X,Y,Z): [%.6f, %.6f, %.6f]\n', ...
    mean(ref_position_skewness(1,:)), mean(ref_position_skewness(2,:)), mean(ref_position_skewness(3,:)));
fprintf('Kurtosis (X,Y,Z): [%.6f, %.6f, %.6f]\n', ...
    mean(ref_position_kurtosis(1,:)), mean(ref_position_kurtosis(2,:)), mean(ref_position_kurtosis(3,:)));
