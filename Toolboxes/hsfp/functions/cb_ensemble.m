 function [cb_ensemble, wts,ens_vals] = cb_ensemble(signal_TT, alpha, prior, target_method)
% Calculates the Conditioned Bayesian Ensemble Posterior (CB_ensemble) 
% flexible probabilities (also refered to as Degree of Conditioning & 
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

% target_method - defines how target value should be selected
% Options: 'latest' - most recent observation of sig at time T (Default).
%          'mean' - average over time T for each sig.         
% (Type: str|char)

%%
% check default method set for z_target
if nargin < 4|| isempty(target_method)
    target_method = 'latest'; 
end

Q = width(signal_TT);

Post_pr_storage = zeros(size(signal_TT));
ens_vals = zeros(1, Q);

% Obtain the posterior Pr & ENS per Signal
for sig = 1:width(signal_TT)
    signal_series = signal_TT{:,sig};
    if strcmp(target_method, 'mean')
        z_target = mean(signal_series);
    elseif strcmp(target_method, 'latest')
        z_target = signal_series(end);
    end
    %% 1. Get posterior via Entropy Pooling
    Post_pr_storage(:,sig) = ep_probs(signal_series, alpha, z_target, prior);
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

%% 6. Weights
weighted_pr_storage = zeros(size(signal_TT));
products = ens_vals.*d_indicators;
wts = products/sum(products);

for sig = 1:Q
weighted_pr_storage(:,sig) = wts(sig)*Post_pr_storage(:,sig);
end

% Sum across the rows
cb_ensemble = sum(weighted_pr_storage, 2)';
end
