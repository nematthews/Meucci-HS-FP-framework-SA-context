function [p_cr, z_ub, z_lb, z_max,z_min] = cr_probs(signal_series, alpha, z_target)
% Single State variable Market conditioning using Crisp Probabilities for
% the HS-FP framework.
% Optimal Prs are set to 1 if z_t falls within the bandwidth range around
% the selected target value.

% INPUT:
% signal_series - smoothed & standardised timeseries for a single state var
% (type: array double, [T x 1])

% alpha - range of probability for bandwidth
% (type: double)

% z_target - target value specific to the signal_series
%          - Options: scalar value chosen by user
%                     'latest'  give latest value in signal_series
%                     'mean' gives average across singal series
% (type: double| str| char)

%%
%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(z_target, 'mean')
        z_target = mean(signal_series);

    elseif strcmp(z_target, 'latest')
        z_target = signal_series(end);
end

%% CDF estimate
try
    Z_sig = signal_series';
    T = length(Z_sig);
    % catch to change possible timetable object to array
catch
    Z_sig = table2array(signal_series)';
    T = length(Z_sig);
end

% Need to get f_Z(z)dz to be able to get leeway range
sorted_z = unique(sort(Z_sig, 'ascend')');   % [1x190 double]


% Sorted state variable cdf
F_est = unique(sum(repmat((1/T), length(sorted_z), 1) .* (Z_sig >= sorted_z), 2));

% Evaluate cdf at target value
cdf_val = interp1(sorted_z, F_est, z_target, 'linear', 'extrap');

%% Quantiles %%%
%%% UPPER
% upper and lower quantiles
z_max = quantile(sorted_z, 1-(alpha/2));
% storage
z_ub = zeros(length(z_target), 1);

%%% LOWER
z_min = quantile(sorted_z, alpha/2);
% storage
z_lb = zeros(length(z_target), 1);

% Pr storage
p_t = zeros(length(z_target), T);
%% Probability Boundary checks
% Check range of values and replace with Max or Min where needed
for i = 1:length(z_target)
    % if cdf @ z*  is lower than alpha/2 bring up to alpha/2
    cdf_val(cdf_val <= alpha/2) = alpha/2;
    % if cdf @ z*  is greater than 1- alpha/2 bring down to 1- alpha/2
    cdf_val(cdf_val >= 1 - alpha/2) = 1 - alpha/2;

    %% Target Boundary checks
    z_star = z_target(i);
    % Adjust for min if out of range
    if z_star < z_min
        z_lb(i) = min(Z_sig);
        z_ub(i) = quantile(Z_sig, cdf_val(i) + alpha/2);
        % Adjust for max if out of range
    elseif z_star > z_max
        z_lb(i) = quantile(Z_sig, cdf_val(i) - alpha/2);
        z_ub(i) = max(Z_sig);

    else
        %% retreive z value @ target_Pr +/- alpha/2
        z_bounds = quantile(Z_sig, [cdf_val(i) - alpha/2, cdf_val(i) + alpha/2]);
        z_lb(i) = z_bounds(1);
        z_ub(i) = z_bounds(2);
    end

    %% Flexible Probability assignment %%%%%
    % Assigning 1 or 0 if z is in or out of R respectively:
    % logical indexing: assign 1 to the elements where the condition is true
    p_t(i, :) = (Z_sig >= z_lb(i) & Z_sig <= z_ub(i));

    % Ensure they sum to 1
    p_cr = zeros(length(z_target), T);
    p_cr(i, :) = p_t(i, :) ./ sum(p_t(i, :));

end


end


