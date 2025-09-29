function [mixture_skewness, mixture_kurtosis] = calculate_mixture_moments(gaussian_means, gaussian_covs, weights, numTimeSteps)
    % Calculate skewness and kurtosis for Gaussian mixture distribution
    % Inputs:
    %   gaussian_means: (3 x numTimeSteps x N_components) - means of each component
    %   gaussian_covs: (3 x 3 x numTimeSteps x N_components) - covariances of each component  
    %   weights: (1 x N_components) - mixture weights
    %   numTimeSteps: number of time steps
    
    N_components = length(weights);
    mixture_skewness = zeros(3, numTimeSteps);
    mixture_kurtosis = zeros(3, numTimeSteps);
    
    for t = 1:numTimeSteps
        % Calculate mixture mean
        mixture_mean = zeros(3, 1);
        for i = 1:N_components
            mixture_mean = mixture_mean + weights(i) * gaussian_means(:, t, i);
        end
        
        % Calculate mixture covariance
        mixture_cov = zeros(3, 3);
        for i = 1:N_components
            diff = gaussian_means(:, t, i) - mixture_mean;
            mixture_cov = mixture_cov + weights(i) * (gaussian_covs(:, :, t, i) + diff * diff');
        end
        
        % Calculate higher-order moments for each dimension
        for dim = 1:3
            sigma2_mix = mixture_cov(dim, dim);
            
            if sigma2_mix < 1e-12  % Avoid numerical issues
                mixture_skewness(dim, t) = 0;
                mixture_kurtosis(dim, t) = 3; % Normal kurtosis for degenerate case
                continue;
            end
            
            % Third central moment (for skewness)
            mu3 = 0;
            for i = 1:N_components
                mu_i = gaussian_means(dim, t, i);
                sigma2_i = gaussian_covs(dim, dim, t, i);
                delta_i = mu_i - mixture_mean(dim);
                
                % For Gaussian component: E[(X_i - μ_mix)^3] = δ^3 + 3δσ_i^2
                % where δ = μ_i - μ_mix
                mu3 = mu3 + weights(i) * (delta_i^3 + 3 * delta_i * sigma2_i);
            end
            
            % Fourth central moment (for kurtosis)  
            mu4 = 0;
            for i = 1:N_components
                mu_i = gaussian_means(dim, t, i);
                sigma2_i = gaussian_covs(dim, dim, t, i);
                delta_i = mu_i - mixture_mean(dim);
                
                % For Gaussian component: E[(X_i - μ_mix)^4] = 
                % δ^4 + 6δ^2σ_i^2 + 3σ_i^4
                % Note: This is different from your original formula
                mu4 = mu4 + weights(i) * (delta_i^4 + 6 * delta_i^2 * sigma2_i + 3 * sigma2_i^2);
            end
            
            % Calculate skewness (standardized third moment)
            mixture_skewness(dim, t) = mu3 / (sigma2_mix^(3/2));
            
            % Calculate kurtosis (standardized fourth moment)
            mixture_kurtosis(dim, t) = mu4 / (sigma2_mix^2);
            
        end
    end
end