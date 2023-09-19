function[srStd,mu3,mu4,T]= sr_std(returns, T, mu3, mu4, sr)
% Calculate the standard deviation of the Sharpe ratio estimator
% distrubution given non-normal returns.
%
% This formula generalizes for both normal and non-normal returns.
% https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1821643
%
%% INPUTS:
%
% 1. returns -  series of differencial/excess realised returns of an
% asset or portfolio
% (type: array double, [T x 1] | timetable object)
%
% 2. T -  Number of returns samples used for `skew`, `kurtosis` and `sr`.
% (Type: scalar)
%
% 3. mu3 - distrubution skewness as third moment expressed in the same
% frequency as the other parameters. NOTE: `skewness`=0 for normal returns.
%
% 4. mu4 - distrubution kurtosis as fourth moment expressed in the same
% frequency as the other parameters. NOTE: `kurtosis`=3 for normal returns.
%
% 4. sr -  Sharpe ratio in the same frequency as the other parameters.

% Author: Nina Matthews (2023)

% $Revision: 1.0 $ $Date: 2023/18/09 16:06:41 $ $Author: Nina Matthews $

%%
% If returns is not a array, convert it to one
if istimetable(returns)
    returns = table2array(returns);
end


% If n is not provided, use the number of samples
if nargin < 2
    T = size(returns, 1);
end

% If skew is not provided, calculate it
if nargin < 3
    mu3 = skewness(returns);
end

% If kurtosis is not provided, calculate it
if nargin < 4
    mu4 = kurtosis(returns,1);
end

% If sr is not provided, estimate it
if nargin < 5
    sr = geo_sr(returns);
end

% % Calculate the standard deviation of the Sharpe ratio
srStd = sqrt((1 + (0.5 * sr.^2) - (mu3 .* sr) + (((mu4 - 3) / 4) * sr.^2)) / (T - 1));


end









