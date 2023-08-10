function [cb_ensemble, wts,ens_vals] = cb_ensemble(signal_TT, alpha, tau_prior, z_target, wt_method)
% Calculates the Conditioned Bayesian Ensemble Posterior (CB_ensemble) 
% flexible probabilities (also referred to as Degree of Conditioning & 
% Correlation (DCC) Flex Prs) via entropy pooling given Q conditioning 
% variables (signals).

%%% NOTE: z_target default for signal q is set as most recent obs of q.
% Code can later be updated to accommodate user defined targets per signal.
% Currently only 'mean' or 'latest' can be specified as method of selecting
% the target.

% INPUTS:
% singal_TT - time series of conditioning variables acting as state signals
% (Type: Timetable, [T x Q])

% alpha - range of probability for bandwidth
% (type: double)

% prior - desired prior probability distribution (usually exp smoothed Prs)
% (type: array double [T x 1])

% z_target - defines how target value should be selected
% Options: 'latest' - most recent observation of sig at time T (Default).
%          'mean' - average over time T for each sig.         
% (Type: str|char)     
% NOTE: unlike other Fns, here a scalar will not be accepted as z_target

% w_method - defines how the final weighted average is calculated 
% Options: 'log-linear' - log-linear weighted average (Default)
%          'simple' - simple weighted average        
% (Type: str|char)

%% check inputs

% check default method set for z_target
if nargin < 4|| isempty(z_target)
    z_target = 'latest'; 
end

% check default weighting method to log-linear
if nargin < 5|| isempty(wt_method)
    wt_method = 'log-linear';
end

%% Entropy & ENS
Q = width(signal_TT);

Post_pr_storage = zeros(size(signal_TT));
ens_vals = zeros(1, Q);

% Obtain the posterior Pr & ENS per Signal
for sig = 1:width(signal_TT)
    signal_series = signal_TT{:,sig};
    if strcmp(z_target, 'mean')
        z_target_val = mean(signal_series);
    elseif strcmp(z_target, 'latest')
        z_target_val = signal_series(end);
    end
    %% 1. Get posterior via Entropy Pooling
    Post_pr_storage(:,sig) = ep_probs(signal_series, alpha, tau_prior, z_target_val);
    %% 2. ENS
    ens_vals(sig) = ens(Post_pr_storage(:,sig)');
end

%% 3. Bhattacharyya coefficient
B_Coeff = zeros(Q, Q);

% Calculate coefficient pairwise across Posterior Prs
for col1 = 1:Q
    for col2 = col1+1:Q
        % Select the columns for pairwise comparison
        column1 = Post_pr_storage(:, col1);
        column2 = Post_pr_storage(:, col2);

        % Calculate the Bhattacharyya coefficient
        bc = sum(sqrt(column1 .* column2));

        % Store the coefficient in the symmetric matrix
        B_Coeff(col1, col2) = bc;
        B_Coeff(col2, col1) = bc;
    end
end

%% 4. Hellinger distance

% Initialize an array to store the Hellinger distances
hellinger_dists = zeros(1, Q);

% Calculate Hellinger distance for each pairwise Bhattacharyya coefficient
for col1 = 1:Q
    for col2 = col1+1:Q
        % Select the columns for pairwise comparison
        bCoeff = B_Coeff(col1, col2);
        
        % Calculate the Hellinger distance using the formula: Hellinger_distance = sqrt(1 - B)
        hellinger_dists(col1, col2) = sqrt(1 - bCoeff);
        hellinger_dists(col2, col1) = sqrt(1 - bCoeff);
    end
end

%% 5. Diversity Indicators
d_indicators = zeros(1, Q);
for sig = 1:Q
    d_indicators(sig) = (1/(Q-1)).*sum(hellinger_dists(:,sig),1);
end

%% 6. Weights using simple weighted average
% weighted_pr_storage = zeros(size(signal_TT));
products = ens_vals.*d_indicators;
wts = products/sum(products);

%% 7. Weighting all Q probabilities together

% 7.1 Simple weighted average
if strcmp(wt_method, 'simple')
    weighted_pr_storage = wts .* Post_pr_storage;
    cb_ensemble = sum(weighted_pr_storage, 2)';
% 7.2 Log-linear weighted average
elseif strcmp(wt_method, 'log-linear')
    logged_priors = log(Post_pr_storage);
    weighted_pr_storage = wts.*logged_priors;
    cb_ensemble = sum(weighted_pr_storage,2)';
    cb_ensemble = exp(cb_ensemble);
end
%%
   % Scale 
    cb_ensemble = cb_ensemble./sum(cb_ensemble);

end


