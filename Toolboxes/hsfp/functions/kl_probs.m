function p_knl = kl_probs(signal_series, h, gamma, z_target)
% Single State variable Market conditioning using Kernel-based Probabilities
% for the HS-FP framework. 
% Optimal Prs are proportional to the distance that z_t is away from the 
% target z^*.

% INPUT:
% signal_series - smoothed & standardised timeseries for a single state var
% (type: array double, [T x 1])

% h - controls the volatility and therefore level of smoothing
% (type: double)

% gamma - determines the tails of the kernel 
% (gamma = 1 -> exp)
% (gamma = 2 -> Gaussian) 
% (type: double)

% z_target - target value specific to the signal_series
% (type: double)
%%
%%%%%%%%%%%%%%%%%%%%%%%
z = signal_series;
 
    % kernel probabilites
    p = exp((-abs(z-z_target).^gamma)/h); 
% Rescale
p_knl = p/sum(p);

end