function pw= ivp(cov)
% ivp performs simple inverse-variance allocation:
% Inverse Variance Portfolio (IVP)

% Copyright 2019 The MathWorks, Inc.

pw= 1./diag(cov);
pw = pw/sum(pw);
end