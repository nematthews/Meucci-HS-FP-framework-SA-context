function EW_Pr_ensemble = ew_ensemble(signal_TT, alpha, prior, target_method)
% Calculates the equally weighted posterior ensemble flexible probabilities
% via entropy pooling given Q conditioning variables (signals).

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
EW_prs_0 = ones(1, Q);
EW_prs = EW_prs_0/sum(EW_prs_0);


weighted_pr_storage = zeros(size(signal_TT));

% Obtain the weighted series of each signal
for sig = 1:width(signal_TT)
    signal_series = signal_TT{:,sig};
    if strcmp(target_method, 'mean')
        z_target = mean(signal_series);
    elseif strcmp(target_method, 'latest')
        z_target = signal_series(end);
    end
    post_pr_T = ep_probs(signal_series, alpha, z_target, prior);
    weighted_pr_storage(:,sig) = EW_prs(sig)*post_pr_T;
end

% Sum across the rows
EW_Pr_ensemble = sum(weighted_pr_storage, 2)';
end






