function [geo_ave,f] = geo_ave(returns,f)
% Calculate the geometric average of input data, allows for annualisation
% or adjusting for other time periods using f. Typically used to calculate
% portfolio returns for annualising, if so "f" = 12 if returns are monthly. 
%
%
%% INPUTS:
% returns: object consisting of J assets and their returns 
% over time period T
% (Type: timetable [T x J])
%
% f - (Optional) number of periods within a yr (e.g 12 for monthly, 
% 4 for quarterly) - Default = 1 i.e not annualised returns
% (type: scalar)  

% Author: Nina Matthews (2023)

% $Revision: 1.1 $ $Date: 2023/013/09 14:18:23 $ $Author: Nina Matthews $

%% Arg specification
if istimetable(returns)
    returns = table2array(returns);
end

if nargin < 2
   f = 1;
end

%% Calculation

% Storage for geometric ave of returns
geo_ave= zeros(1,width(returns));

for i = 1:width(returns)

    PRet = returns(:,i);
    % Number of periods under analysis
    n = height(PRet);
    % Calculate the geometric mean return (in decimals)
    geo_ave(i) = ((prod(1 + PRet))^ (f / n)) - 1;
end

end

