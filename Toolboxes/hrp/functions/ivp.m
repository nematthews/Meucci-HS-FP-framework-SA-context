function pw= ivp(AssetCov)
% ivp performs simple inverse-variance allocation:
% Inverse Variance Portfolio (IVP)
%INPUTS:
% Assectov - covaraince matrix of asset returns


%%
pw= 1./diag(AssetCov);
pw = pw/sum(pw);
end