function expected_max_sr = est_max_sr(independent_trials, trials_sr_Var, expected_mean_sr)
% Computes the expected maximum Sharpe ratio analytically, typically used
% to estimate the DSR. 
%
%% INPUT:
%
%  trials_returns -  Returns from all trials.
% (type: array double, [T x M] | timetable object)
%
%  expected_mean_sr -  expected mean Sharpe ratio (in same frequency 
% as the psr, typically 0 for random strategies.
% (type: double)
% 
%  independent_trials: Number of independent trials E[K] estimated from 
%  kmeans clustering (see clusterKMeansBase, clusterKMeansTop)
% (type: double, between 1 and E[K]).
% 
%  trials_sr_std: Standard deviation of estimated Sharpe ratios for all trials.
% (type: double)


% $Revision: 1.2 $ $Date: 2023/09/20 10:46:01 $ $Author: Nina Matthews $

%% 
  
   % Input checks:

    % if nargin < 2
    %     independent_trials = num_independent_trials(trials_returns);
    % end
    % 
    % if nargin < 3
    %     srs = estimated_sharpe_ratio(trials_returns);
    %     trials_sr_std = std(srs);
    % end

    if nargin < 3
        expected_mean_sr = 0;
   end

    % Define the Euler-Mascheroni constant
    emc = 0.5772156649; 

    maxZ = (1 - emc) * norminv(1 - 1/independent_trials) + emc * norminv(1 - 1/(independent_trials * exp(1)));
    expected_max_sr = expected_mean_sr + (sqrt(trials_sr_Var) * maxZ);
end

% function num_trials = num_independent_trials(trials_returns)
%     % Your implementation for num_independent_trials function goes here
%     % Implement based on your requirements
% end
% 
% function srs = estimated_sharpe_ratio(trials_returns)
%     % Your implementation for estimated_sharpe_ratio function goes here
%     % Implement based on your requirements
% end
