function p_rw = rw_probs(w_len)

% rw_probs facilitates the calculation of the Time-conditioned probabilities
% used in the HS-FP framework developed by Meucci (2010).
% It uses window lengths specified in number of months to generate
% the sequence of normalised Prs. Default window length = 60 mnths (5 yrs).

% check default value: 60 = 5 yr rolling window
if nargin < 1 || isempty(w_len)
    w_len = 60; % Default value of w_len is set to 60
end

% p(end-w_len:end) = 1;
% sets first w_len elements in the p zero array of length T = to 1
p(1:w_len) = 1;
% normalises the array p so that the sum of its elements = 1.
p_rw = p / sum(p);
end

