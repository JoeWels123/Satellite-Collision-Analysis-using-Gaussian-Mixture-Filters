function [log_collision_probs, collision_probs, ratios, sigma_x_ratios, sigma_z_ratios] = compute_collision_probability_over_time_weighted(...
    means1_pos, covs1_pos, vels1, ...
    means2_pos, covs2_pos, vels2, R)

T = size(means1_pos, 2); % Number of time steps
log_collision_probs = zeros(1, T);
collision_probs = zeros(1, T);
ratios = zeros(1, T);
sigma_x_ratios = zeros(1, T);
sigma_z_ratios = zeros(1, T);

for i = 1:T
    % Extract states and covariances at time step i
    mu1 = means1_pos(:, i); % [3x1] ISS position
    mu2 = means2_pos(:, i); % [3x1] Debris position
    P1 = covs1_pos(:, :, i); % [3x3] ISS position covariance
    P2 = covs2_pos(:, :, i); % [3x3] Debris position covariance
    v1 = vels1(:, i); % [3x1] ISS velocity
    v2 = vels2(:, i); % [3x1] Debris velocity

    % Relative quantities
    rel_mu = mu1 - mu2;
    rel_vel = v1 - v2;
    
    if norm(rel_vel) < 1e-8
        % Objects have nearly identical velocities - high collision risk
        log_collision_probs(i) = 0; % log(1) = 0, meaning Pc = 1
        collision_probs(i) = 1;
        ratios(i) = 0;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    % Construct encounter frame
    ye = rel_vel / norm(rel_vel);
    x_temp = rel_mu - dot(rel_mu, ye) * ye;
    
    if norm(x_temp) < 1e-12
        % Objects are on collision course
        log_collision_probs(i) = 0; % log(1) = 0, meaning Pc = 1
        collision_probs(i) = 1;
        ratios(i) = 0;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    xe = x_temp / norm(x_temp);
    ze = cross(xe, ye);
    Me = [xe, ye, ze]; % 3x3 encounter frame transformation matrix
    
    % Transform covariances to encounter frame
    Pe1 = Me' * P1 * Me;
    Pe2 = Me' * P2 * Me;
    Pe = Pe1 + Pe2;
    
    % Drop along-track (ye) component → 2D in encounter plane
    Pe_2D = Pe([1 3], [1 3]); % xe-ze plane
    mu_rel_2D = Me' * rel_mu;
    mu_2D = mu_rel_2D([1 3]);
    
    % Check if covariance matrix is positive definite
    if det(Pe_2D) <= 0 || any(eig(Pe_2D) <= 0)
        log_collision_probs(i) = -Inf;
        collision_probs(i) = 0;
        ratios(i) = Inf;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    % Rotate covariance to principal axes
    [V, D] = eig(Pe_2D);
    P_diag = D;
    mu_rot = V' * mu_2D;
    sigma_x = sqrt(abs(P_diag(1,1)));
    sigma_z = sqrt(abs(P_diag(2,2)));
    dx = mu_rot(1);
    dz = mu_rot(2);
    
    % Check for numerical issues
    if sigma_x < 1e-12 || sigma_z < 1e-12
        log_collision_probs(i) = -Inf;
        collision_probs(i) = 0;
        ratios(i) = Inf;
        sigma_x_ratios(i) = NaN;
        sigma_z_ratios(i) = NaN;
        continue;
    end
    
    % Compute collision probability using double integral
    % Define the bivariate normal PDF in rotated coordinates
    pdf_func = @(x, z) (1 / (2 * pi * sigma_x * sigma_z)) * ...
               exp(-0.5 * ((x - dx).^2 / sigma_x^2 + (z - dz).^2 / sigma_z^2));
    
    % Integration limits for circular region
    z_lower = @(x) -sqrt(max(0, R^2 - x.^2));
    z_upper = @(x) sqrt(max(0, R^2 - x.^2));
    
    try
        % Compute the double integral
        Pc = integral2(pdf_func, -R, R, z_lower, z_upper, ...
                      'AbsTol', 1e-12, 'RelTol', 1e-9);
        
        % Handle numerical precision issues
        Pc = max(0, min(1, Pc)); % Clamp to [0, 1]
        
        if Pc > 0
            log_Pc = log(Pc);
        else
            log_Pc = -Inf;
        end
        
    catch ME
        % If integration fails, fall back to very small probability
        warning('Integration failed at time step %d: %s', i, ME.message);
        Pc = 0;
        log_Pc = -Inf;
    end
    
    % Compute ratio to check approximation validity
    miss_distance = norm(mu_2D);
    uncertainty_radius = sqrt(sigma_x^2 + sigma_z^2);
    ratio = miss_distance / uncertainty_radius;
    
    % Compute sigma/R ratios
    sigma_x_ratios(i) = sigma_x / R;
    sigma_z_ratios(i) = sigma_z / R;
    
    % Store results
    log_collision_probs(i) = log_Pc;
    collision_probs(i) = Pc;
    ratios(i) = ratio;
end
end