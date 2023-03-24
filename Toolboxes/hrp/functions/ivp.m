function pw= ivp(AssetCov)
% ivp performs simple inverse-variance allocation:
% Inverse Variance Portfolio (IVP)

% Copyright 2019 The MathWorks, Inc.

pw= 1./diag(AssetCov);
pw = pw/sum(pw);
end