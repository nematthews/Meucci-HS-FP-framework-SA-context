function ivp = clusterIVP(cov, use_extended_terms)
    % Check if use_extended_terms is provided, otherwise set it to false
    if nargin < 2
        use_extended_terms = false;
    end
    
    % Compute the minimum-variance portfolio
    ivp = 1 ./ diag(cov);
    
    if use_extended_terms
        n = size(cov, 1);
        corr = cov2corr(cov);
        
        % Obtain average off-diagonal correlation
        rho = (sum(corr(:)) - n) / (n^2 - n);
        
        invSigma = sqrt(ivp);
        ivp = ivp - rho * invSigma * sum(invSigma) / (1 + (n - 1) * rho);
    end
    
    ivp = ivp / sum(ivp);
end
