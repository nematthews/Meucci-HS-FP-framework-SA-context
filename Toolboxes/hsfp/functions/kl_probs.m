function p_knl = kl_probs(signal_series, h, gamma, z_target)
% Single State variable Market conditioning using Kernel-based Probabilities
% for the HS-FP framework.
% Optimal Prs are proportional to the distance that z_t is away from the
% target z^*.
%
%% INPUT:
% signal_series - smoothed & standardised timeseries for a single state var
% (type: array double, [T x 1])
%
% h - controls the volatility and therefore level of smoothing
% (type: double)
%
% gamma - determines the tails of the kernel
% (gamma = 1 -> exp)
% (gamma = 2 -> Gaussian)
% (type: double)
%
% z_target - target value specific to the signal_series
%          - Options: scalar value chosen by user
%                     'latest'  give latest value in signal_series
%                     'mean' gives average across singal series
% (type: double| str| char)


% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/02/20 16:10:46 $ $Author: Nina Matthews $

%%
%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(z_target, 'mean')
    z_target = mean(signal_series);

elseif strcmp(z_target, 'latest')
    z_target = signal_series(end);
end

z = signal_series;

% kernel probabilites
try
    p = exp((-abs(z-z_target).^gamma)/h);
    % catch to change possible timetable object to array
catch
    z = table2array(signal_series);
    p = exp((-abs(z-z_target).^gamma)/h);
end
% Rescale
p_knl = (p/sum(p))';

end