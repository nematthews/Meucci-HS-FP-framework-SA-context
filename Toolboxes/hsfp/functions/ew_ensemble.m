function [EW_Pr_ensemble] = ew_ensemble(signal_TT, alpha, tau_prior, z_target)
% Calculates the equally weighted posterior ensemble flexible probabilities
% via entropy pooling given Q conditioning variables (signals).
%
%%% NOTE: z_target default for signal q is set as most recent obs of q.
% Code can later be updated to accommodate user defined targets per signal.
% Currently only 'mean' or 'latest' can be specified as method of selecting
% the target.
%
%% INPUTS:
% singal_TT - time series of conditioning variables acting as state signals
% (Type: Timetable, [T x Q])
%
% alpha - range of probability for bandwidth
% (type: double)
%
% tau_prior - tau of desired prior probability distribution (exp smoothed Prs)
% (type: double)
%
% z_target - defines how target value should be selected
% Options: 'latest' - most recent observation of sig at time T (Default).
%          'mean' - average over time T for each sig.         
% (Type: str|char)     
% NOTE: unlike other Fns, here a scalar will not be accepted as z_target


% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/02/20 16:10:46 $ $Author: Nina Matthews $

%% Input checks
% check default method set for z_target
if nargin < 4|| isempty(z_target)
    z_target = 'latest'; 
end
%%
Q = width(signal_TT);
EW_prs_0 = ones(1, Q);
EW_prs = EW_prs_0/sum(EW_prs_0);


weighted_pr_storage = zeros(size(signal_TT));

% Obtain the weighted series of each signal
for sig = 1:width(signal_TT)
    signal_series = signal_TT{:,sig};
    if strcmp(z_target, 'mean')
        z_target_val = mean(signal_series);
    elseif strcmp(z_target, 'latest')
        z_target_val = signal_series(end);
    end
    post_pr_T = ep_probs(signal_series, alpha, tau_prior, z_target_val);
    weighted_pr_storage(:,sig) = EW_prs(sig)*post_pr_T;
end

% Sum across the rows
EW_Pr_ensemble = sum(weighted_pr_storage, 2)';
end






