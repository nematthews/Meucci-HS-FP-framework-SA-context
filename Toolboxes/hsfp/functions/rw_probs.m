function [p_rw,w_len] = rw_probs(Sig_timeTable,w_len)

% rw_probs facilitates the calculation of the Time-conditioned probabilities
% used in the HS-FP framework developed by Meucci (2010).
% It uses window lengths specified in number of months to generate
% the sequence of normalised Prs. Default window length = 60 mnths (5 yrs).
%
%%INPUTS
% Sig_timeTable - time series of state signals (any frequency)
% (Type: Timetable object | array)
%
% w_len - window length (in same frequency as X_timeTable)
% (Type: double)


% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/02/20 16:10:46 $ $Author: Nina Matthews $

%%
% check default value: 60 = 5 yr rolling window
if nargin < 2 || isempty(w_len)
    w_len = 60; % Default value of w_len is set to 60
end

% p(end-w_len:end) = 1;
% sets first w_len elements in the p zero array of length T = to 1
p(1:w_len) = 1;
% normalises the array p so that the sum of its elements = 1.
p_rw = p / sum(p);

if nargin > 0 && ~isempty(Sig_timeTable)
    % Adjust the length of p_rw to include 0's outside of window
    x = length(p_rw);
    y = height(Sig_timeTable);
    
    if x < y
        zerosToAdd = y - x;
        p_rw = [zeros(1, zerosToAdd), p_rw];
    end
end 
end

