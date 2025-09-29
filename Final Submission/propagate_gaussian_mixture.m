function [collisionprobs_total, lognoprob, noprob] = propagate_gaussian_mixture(...
    Object1_gaussian_means, Object2_gaussian_means, ...
    Object1_gaussian_covariances, Object2_gaussian_covariances, ...
    Object1_gaussian_weights, Object2_gaussian_weights, ...
    numTimeSteps, N, R)

global common_timeVector

% Initialize arrays to store collision probabilities
collisionprobs_total = zeros(1, numTimeSteps);

% Initialize total collision probability arrays
lognoprob = zeros(1, numTimeSteps);
noprob = zeros(1, numTimeSteps);

% Initialize collision probability for this time series
collision_prob_timestep = zeros(1, numTimeSteps);

% Loop over all combinations of Gaussian components
for g1 = 1:N  % Object1 Gaussians
    for g2 = 1:N  % Object2 Gaussians
        
        % Extract means and covariances for this Gaussian pair
        ISS_gaussian_mean = Object1_gaussian_means(:, :, g1);  % [6 x numTimeSteps]
        debris_gaussian_mean = Object2_gaussian_means(:, :, g2);  % [6 x numTimeSteps]
        ISS_gaussian_cov = Object1_gaussian_covariances(:, :, :, g1);  % [6 x 6 x numTimeSteps]
        debris_gaussian_cov = Object2_gaussian_covariances(:, :, :, g2);  % [6 x 6 x numTimeSteps]
        
        % Extract position and velocity components
        ISS_pos = ISS_gaussian_mean(1:3, :);     % [3 x numTimeSteps]
        ISS_vel = ISS_gaussian_mean(4:6, :);     % [3 x numTimeSteps]
        debris_pos = debris_gaussian_mean(1:3, :);  % [3 x numTimeSteps]
        debris_vel = debris_gaussian_mean(4:6, :);  % [3 x numTimeSteps]
        
        ISS_pos_cov = ISS_gaussian_cov(1:3, 1:3, :);     % [3 x 3 x numTimeSteps]
        debris_pos_cov = debris_gaussian_cov(1:3, 1:3, :);  % [3 x 3 x numTimeSteps]
        
        % Compute collision probability for this Gaussian pair
        [~, collisionprobs_pair, ratios] = compute_collision_probability_over_time_weighted(...
            ISS_pos, ISS_pos_cov, ISS_vel, ...
            debris_pos, debris_pos_cov, debris_vel, R);
        
        % Weight by the product of Gaussian weights and add to total
        weight_product = Object1_gaussian_weights(g1) * Object2_gaussian_weights(g2);
        collision_prob_timestep = collision_prob_timestep + weight_product * collisionprobs_pair;
        
        % Optional: Print progress for debugging
        if mod(g1*N + g2, 10) == 0
            fprintf('Processed Gaussian pair (%d,%d), max Pc for this pair: %.6e\n', ...
                g1, g2, max(collisionprobs_pair));
        end
    end
end

% Store the total collision probability
collisionprobs_total(1, :) = collision_prob_timestep;

% Create validity check plot for the mixture result
figure('Position', [100, 100, 1200, 600]);
subplot(2,1,1);
semilogy(common_timeVector, collision_prob_timestep, 'b-', 'LineWidth', 2);
hold on;

high_prob_threshold = 1e-4;
high_prob_indices = collision_prob_timestep > high_prob_threshold;

if any(high_prob_indices)
    semilogy(common_timeVector(high_prob_indices), collision_prob_timestep(high_prob_indices), ...
             'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
end

yline(high_prob_threshold, 'm--', 'LineWidth', 1.5, 'Label', 'Pc = 10^{-4}');
xlabel('Time');
ylabel('Collision Probability');
title('Mixture-Based Collision Probability over Time');
grid on;

subplot(2,1,2);
% For the mixture, we can't easily compute a single "ratio" metric
% Instead, show the collision probability on linear scale
plot(common_timeVector, collision_prob_timestep, 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('Collision Probability (Linear Scale)');
title('Collision Probability - Linear Scale');
grid on;

fprintf('Max collision probability = %.6e\n', max(collision_prob_timestep));

end

