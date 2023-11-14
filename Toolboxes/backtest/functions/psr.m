function psr = psr(dif_returns, sr_benchmark, sr, srStd)
% Calculate the Probabilistic Sharpe Ratio (PSR) - Bailey & de Prado (2012)
% PSR(SR*) = probability that SR^ > SR* where, SR^ = estimated sharpe
% ratio, and SR* is a benchmark comparative strategy.
%
%% INPUT:
%
% dif_returns - series of differential/excess realised returns of an 
% asset or portfolio.
% (type: array double, [T x 1] | timetable object)
%
% NOTE: dif_returns can be a series of an asset returns, in this case the
% excess return used benchmark return as 0 as measure of comparing against
% no investment skill.
%
% sr_benchmark - Benchmark SR expressed in the same frequency as other parameters.
% By default set to zero to compare strategy to "no investment skill'. 
% Vector can be supplied of benchmark strategy to calc sr_benchmark.
% (type: scalar | double [T x 1] | timetable object)
%
% sr - Sharpe ratio in the same frequency as other parameters. (calculated 
% if not provided)
% (type: scalar)  
%
% srStd -  Standard deviation of the Estimated SR. (calculated if not
% provided)
% (type: scalar)  
%
% NOTE:
% -----
% PSR(SR*) = probability that SR^ > SR*
% SR^ = sharpe ratio estimated with "dif_returns", or "sr"
% SR* = "sr_benchmark"

% Author: Nina Matthews (2023)


% $Revision: 1.2 $ $Date: 2023/09/25 14:03:01 $ $Author: Nina Matthews $
%% Input checks
if istimetable(dif_returns)
    dif_returns = table2array(dif_returns);
end

if istimetable(sr_benchmark)
    sr_benchmark = table2array(sr_benchmark);
end


% This is the case where benchmark = 0 as measure of comparing against
% no investment skill. (IF sr or srStd is not provided)
if nargin < 2
    sr_benchmark = 0;
end

% check if sr_benchmark is a scalar or if geo_sr needs to be calculated
% from array
[r,c] = size(sr_benchmark);
if (r > 1) || (c > 1)
   sr_benchmark = geo_sr(sr_benchmark,1); 
end

% If sr and srStd are not provided, estimate them from returns
if nargin < 3
    % Specify f = 1 to ensure not to annualise
    sr = geo_sr(dif_returns,1);
end
disp(sr)

if nargin < 4
    srStd = sr_std(dif_returns);
end

%% Calculate the PSR using the normal cumulative distribution function
psr = normcdf((sr - sr_benchmark) / srStd);


end



