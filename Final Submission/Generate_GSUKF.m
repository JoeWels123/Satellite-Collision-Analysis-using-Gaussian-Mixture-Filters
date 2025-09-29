function [ISS_mixture_cov, ISS_mixture_mean, ISS_gaussian_means, ISS_gaussian_covs, ISS_mu_components, ISS_P_components, gaussian_weights, Wm_all, Wc_all] = Generate_GSUKF(mu, P, N, numSamples, startTime, stopTime, sampleTime, File)
    timeVector = startTime:seconds(sampleTime):stopTime;
    numTimeSteps = length(timeVector);

    % Unscented transform parameters
    alpha = 1;
    beta = 2;
    kappa = -3;
    n = 6;
    n_sigma = 2*n + 1;
    
    % Generate samples and fit Gaussian Mixture Model
    data = mvnrnd(mu, P, numSamples);


    options = statset('MaxIter', 10000, 'Display', 'off');
    fprintf('Number of samples: %d\n', size(data', 1));
    fprintf('Number of dimensions: %d\n', size(data', 2));
    fprintf('Requested GMM components: %d\n', N);
    gm = fitgmdist(data, N, 'Options', options, 'RegularizationValue', 1e-30);
    

    N_components = gm.NumComponents;
    
    % Extract GMM parameters
    ISS_mu_components = gm.mu';
    ISS_P_components = gm.Sigma;
   

    % Set equal weights for all components
    gaussian_weights = gm.ComponentProportion;
    % Initialize storage for sigma points and weights for each component
    ISS_all_sigma_points = zeros(n, n_sigma, N_components);
    Wm_all = zeros(n_sigma, N_components);
    Wc_all = zeros(n_sigma, N_components);
    
    % Generate sigma points for each Gaussian component
    for i = 1:N_components
        mu_i = ISS_mu_components(:, i);
        P_i = ISS_P_components(:, :, i);
        
        [sigma_points_i, Wm_i, Wc_i] = generate_sigma_points(mu_i, P_i, alpha, beta, kappa);
        
        
        keplerian_sigma_points_i = zeros(n, n_sigma);
        for j = 1:n_sigma
            classical_elements = equinoctial_to_classical(sigma_points_i(:, j));
           
            % Extract fields from struct
            keplerian_sigma_points_i(:, j) = [
                classical_elements.meanmotion;      % mean motion
                classical_elements.eccentricity;      % eccentricity
                classical_elements.inclination;      % inclination
                classical_elements.raan;   % right ascension of ascending node
                classical_elements.argp;      % argument of periapsis
                classical_elements.meanAnomaly      % mean anomaly
            ];
           
        end
        

        ISS_all_sigma_points(:, :, i) = keplerian_sigma_points_i;
        Wm_all(:, i) = Wm_i;
        Wc_all(:, i) = Wc_i;
    end
    
    sc_ISS = satelliteScenario(startTime, stopTime, sampleTime);

    ISS_positions = cell(N_components, n_sigma);
    ISS_velocities = cell(N_components, n_sigma);
    disp(startTime)
    for j = 1:N_components
        for i = 1:n_sigma
            elems = struct( ...
                'meanMotion', ISS_all_sigma_points(1, i, j) * 240, ... % rev/day
                'eccentricity', ISS_all_sigma_points(2, i, j), ...
                'inclination', ISS_all_sigma_points(3, i, j), ...
                'raan', ISS_all_sigma_points(4, i, j), ...
                'argPeriapsis', ISS_all_sigma_points(5, i, j), ...
                'meanAnomaly', ISS_all_sigma_points(6, i, j) ...
            );
            
            % Create unique filename for each sigma point
            outputFile = 'updated.xml';
            updateOMMFile(File, outputFile, elems);
            
            sat = satellite(sc_ISS, outputFile);

            [pos, vel, ~] = states(sat);
  
            ISS_positions{j, i} = pos;
            ISS_velocities{j, i} = vel;

            delete(outputFile)
        end
    end
    
    ISS_predicted_mean_cartesian = zeros(6, numTimeSteps, N_components);
    ISS_predicted_cov_cartesian = zeros(6, 6, numTimeSteps, N_components);
    ISS_sigma_points_cartesian = zeros(6, n_sigma, numTimeSteps, N_components);
    
    for j = 1:N_components
        for t = 1:numTimeSteps
            for i = 1:n_sigma
                r_t = ISS_positions{j,i}(:,t);  % This will work after fixing storage
                v_t = ISS_velocities{j,i}(:,t);
                ISS_sigma_points_cartesian(:,i, t, j) = [r_t; v_t];
            end
            
            Wm_j = Wm_all(:,j);
            Wc_j = Wc_all(:,j);
            mean_jt = ISS_sigma_points_cartesian(:, :, t, j) * Wm_j;
            
            cov_jt = zeros(6, 6);
            for i = 1:n_sigma
                diff = ISS_sigma_points_cartesian(:,i, t, j) - mean_jt;
                cov_jt = cov_jt + Wc_j(i) * (diff * diff');
            end
            
            ISS_predicted_mean_cartesian(:,t,j) = mean_jt;
            ISS_predicted_cov_cartesian(:,:,t,j) = cov_jt;
        end
    end
    ISS_predicted_cov_keplerian = zeros(6, 6, N_components);
    
    for j = 1:N_components
        Wm_j = Wm_all(:,j);
        Wc_j = Wc_all(:,j);
        
        mean_jt = ISS_all_sigma_points(:, :, j) * Wm_j; 
        
        cov_jt = zeros(6, 6);
        for i = 1:n_sigma
            diff = ISS_all_sigma_points(:,i, j) - mean_jt;
            cov_jt = cov_jt + Wc_j(i) * (diff * diff');
        end
        
        ISS_predicted_cov_keplerian(:,:,j) = cov_jt;
    end
    ISS_gaussian_means = ISS_predicted_mean_cartesian; 
    ISS_gaussian_covs = ISS_predicted_cov_cartesian;               
    
    ISS_mixture_mean = zeros(6, numTimeSteps);
    ISS_mixture_cov = zeros(6, 6, numTimeSteps);
    
    for t = 1:numTimeSteps
        mean_mix_t = zeros(6,1);
        for j = 1:N_components
            mean_j = ISS_gaussian_means(:,t,j);
            w_j = gaussian_weights(j);
            mean_mix_t = mean_mix_t + w_j * mean_j;
        end
        ISS_mixture_mean(:,t) = mean_mix_t;
        
        cov_mix_t = zeros(6, 6);
        for j = 1:N_components
            mean_j = ISS_gaussian_means(:,t,j); 
            cov_j = ISS_gaussian_covs(:,:,t,j);
            w_j = gaussian_weights(j);
            
            diff = mean_j - mean_mix_t;
    
            cov_mix_t = cov_mix_t + w_j * (cov_j + diff * diff'); 
        end
        ISS_mixture_cov(:,:,t) = cov_mix_t;
    end
end

