function [geo_SR] = geo_sr(dif_returns, f)
% Calculate the Sharpe Ratio based off geometric averages of the risky
% asset and a risk free asset given realised returns of the risky and risk
% less.
%
%% INPUT:
% dif_returns - series of differential/excess realised returns of an asset or portfolio
% (type: array double, [T x 1] | timetable object)
%
% NOTE: dif_returns can be a series of an asset returns, in this case the
% excess return used benchmark return as 0 as measure of comparing against
% no investment skill.
%
% f - (Optional) number of periods within a yr (e.g 12 for monthly, 
% 4 for quarterly) only used to calculate geo_ave (i.e annualises mean
% returns NOT the SR.
% (type: double) NOTE: default = 1 i.e not annualised returns
% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/09/20 10:46:01 $ $Author: Nina Matthews $


if istimetable(dif_returns)
    dif_returns = table2array(dif_returns);
end

% Number of periods within a year (e.g., 12 for monthly) DEFUALT HERE f = 1
if nargin < 2
    f = 1;
end

%% Geometric differential/excess returns for SR
% We use geometric average here as we are not using SR for prediction and
% decision making. We using historical realised excess returns to
% calculate an indicative measure of historical performance.

geo_SR = zeros(1,width(dif_returns));
ExcessRet_sd = zeros(1,width(dif_returns));

for i = 1:width(dif_returns)
    % Calculate Geometric Average of Excess Returns (to annualise f = 12)
    geo_ExcessRet_ave = geo_ave(dif_returns(:,i),f);
    % Calc sd (under the assumption of normality)
    ExcessRet_sd(i) = std(dif_returns(:,i),1);

    %% Calculate SRatio
    geo_SR(i) = geo_ExcessRet_ave/ExcessRet_sd(i);
end
end



