function [z_ub, z_lb, p_cr] = cr_probs(signal_series, alpha, z_target)
% Single State variable Market conditioning using Crisp Probabilities for 
% the HS-FP framework. 
% Optimal Prs are set to 1 if z_t falls within the bandwidth range around 
% the selected target value.

% INPUT:
% signal_series - smoothed & standardised timeseries for a single state var
% (type: array double, [1 x T])

% alpha - range of probability for bandwidth
% (type: double)

% z_target - target value specific to the signal_series
% (type: double)
%%
%%%%%%%%%%%%%%%%%%%%%%%
%% CDF estimate
Z_sig = signal_series';
T = length(Z_sig);

% Need to get f_Z(z)dz to be able to get leeway range
sorted_z = unique(sort(Z_sig, 'ascend')');   % [1x190 double]

% Potential alternative for F:
F_est = unique(sum(repmat((1/T), length(sorted_z), 1) .* (Z_sig >= sorted_z), 2));

% Interpolate to get full est cdf
cdf = interp1(sorted_z, F_est, z_target, 'linear', 'extrap');

%% Quantiles %%%
%%% UPPER
% upper and lower quantiles of z
z_max = quantile(sorted_z, 1 - (alpha/2));
% storage
z_ub = zeros(length(z_target), 1);

%%% LOWER
z_min = quantile(sorted_z, alpha/2);
% storage
z_lb = zeros(length(z_target), 1);

% Pr storage
p_t = zeros(length(z_target), T);
%% Boundary checks
% Check range of values and replace with Max or Min where needed
for i = 1:length(z_target)
    cdf(cdf <= alpha/2) = alpha/2;
    cdf(cdf >= 1 - alpha/2) = 1 - alpha/2;

    sorted_z = z_target(i);
    % Adjust for min
    if sorted_z < z_min
        z_lb(i) = min(Z_sig);
        z_ub(i) = quantile(Z_sig, cdf(i) + alpha/2);
        % Adjust for max
    elseif sorted_z > z_max
        z_lb(i) = quantile(Z_sig, cdf(i) - alpha/2);
        z_ub(i) = max(Z_sig);

    else
        z_bounds = quantile(Z_sig, [cdf(i) - alpha/2, cdf(i) + alpha/2]);
        z_lb(i) = z_bounds(1);
        z_ub(i) = z_bounds(2);
    end

    %% Probability assignment %%%%%
    % Assigning 1 or 0 if z is in or out of R respectively:
    % logical indexing: assign 1 to the elements where the condition is true
    p_t(i, :) = (Z_sig >= z_lb(i) & Z_sig <= z_ub(i));

    % Ensure they sum to 1
    p_cr = zeros(length(z_target), T);
    p_cr(i, :) = p_t(i, :) ./ sum(p_t(i, :));
end


end


