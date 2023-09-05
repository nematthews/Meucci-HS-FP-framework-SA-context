function pw= ivp(AssetCov)
% ivp performs simple inverse-variance allocation:
% Inverse Variance Portfolio (IVP)
%INPUTS:
% Assectov - covaraince matrix of asset returns


% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/02/20 16:10:46 $ $Author: Nina Matthews $


%%
pw= 1./diag(AssetCov);
pw = pw/sum(pw);
end