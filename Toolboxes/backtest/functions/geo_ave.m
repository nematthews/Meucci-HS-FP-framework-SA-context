function [geo_ave,f] = geo_ave(Returns,f)
% Calculate the geometric average of input data, allows for annualisation
% or adjusting for other time periods using f.
%
%
%% INPUTS:
% AssetReturns: object consisting of J assets and their returns 
% over time period T
% (Type: timetable [T x J])
%
% f - Number of periods within a yr (e.g 12 for monthly, 4 for quarterly)
% (type: double) -Defualt  =  1
% Author: Nina Matthews (2023)


% Author: Nina Matthews (2023)

% $Revision: 1.1 $ $Date: 2023/013/09 14:18:23 $ $Author: Nina Matthews $

%% Arg specification
if istimetable(Returns)
    Returns = table2array(Returns);
end

if nargin < 2
   f = 1;
end

%% Calculation

% Storage for geometric ave of returns
geo_ave= zeros(1,width(Returns));

for i = 1:width(Returns)

    PRet = Returns(:,i);
    % Number of periods under analysis
    n = height(PRet);
    % Calculate the geometric mean return (in decimals)
    geo_ave(i) = ((prod(1 + PRet))^ (f / n)) - 1;
end

end

