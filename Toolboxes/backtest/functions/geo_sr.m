function [geo_SR] = geo_sr(SR_dif_t, f)
% Calculate the Sharpe Ratio based off geometric averages of the risky
% asset and a risk free asset given realised returns of the risky and risk
% less.
%
%% INPUT:
% SR_dif_t - series of differencial/excess realised returns of an asset or portfolio
% (type: array double, [T x 1] | timetable object)
%
% rfr_series - series of realised returns on a risk free asset eg cash
% (type: array double, [T x 1])
%
% f - Number of periods within a yr (e.g 12 for monthly, 4 for quarterly)
% (type: double)
% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/09/12 15:03:01 $ $Author: Nina Matthews $


if istimetable(SR_dif_t)
    SR_dif_t = table2array(SR_dif_t);
end

% Number of periods within a year (e.g., 12 for monthly)
if strcmp(f, 'annualise')
    f = 12;
else 
    f = 1;
end

%% Annualised geometric differencial/excess returns for SR

geo_ExcessRet_ave = geo_ave(SR_dif_t,f);
%% Calculate Sigma of differencial/excess returns
SR_dif_sd = std(geo_ExcessRet_ave,1);

%% Calcuate SRatio

geo_SR = geo_ExcessRet_ave/SR_dif_sd;



 

