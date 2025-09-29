function d = bhattacharyya_distance(mu1, Sigma1, mu2, Sigma2)
% Computes the Bhattacharyya distance between two multivariate Gaussians.
% 
% Inputs:
%   mu1, mu2     - Mean vectors (nx1)
%   Sigma1, Sigma2 - Covariance matrices (nxn)
%
% Output:
%   d - Bhattacharyya distance (scalar)

    % Check dimension consistency
    assert(all(size(mu1) == size(mu2)), 'Mean vectors must be the same size.');
    assert(all(size(Sigma1) == size(Sigma2)), 'Covariance matrices must be the same size.');
    
    % Average covariance
    Sigma = 0.5 * (Sigma1 + Sigma2);
    
    % Regularize in case of near-singular matrix
    epsilon = 1e-6;
    Sigma = Sigma + epsilon * eye(size(Sigma));
    
    % First term (Mahalanobis-like term)
    diff = mu1 - mu2;
    term1 = 0.125 * (diff') * (Sigma \ diff);
    
    % Second term (determinant term)
    det_Sigma1 = det(Sigma1 + epsilon * eye(size(Sigma1)));
    det_Sigma2 = det(Sigma2 + epsilon * eye(size(Sigma2)));
    det_Sigma  = det(Sigma);
    
    % Handle numerical issues with determinants
    if det_Sigma <= 0 || det_Sigma1 <= 0 || det_Sigma2 <= 0
        term2 = 0;
    else
        term2 = 0.5 * log(det_Sigma / sqrt(det_Sigma1 * det_Sigma2));
    end
    
    % Bhattacharyya distance
    d = term1 + term2;
end
