function [SR_opt_wts,MV_SR_sdv, MV_SR_return,frontier_rsk_vec, frontier_return_vec] = maxsr(AssetList,AssetMean,AssetCovar, Rfr)
% Using MV asset allocation

% Mean-Variance Opimation with a Maximum Sharpe Ratio objective function
% for portfolio weights allocation.

% INPUTS:
% AssetReturns: object consisting of J assets and their returns 
% over time period T
% (Type: timetable [T x J])

%% Set up data for MV SR portfolio object


%% Set up portfolio object
p = Portfolio('AssetList', AssetList, 'RiskFreeRate', Rfr);
p = setAssetMoments(p, AssetMean, AssetCovar);
% Set initial weights as 1/N portfolio
% p = setInitPort(p, 1/p.NumAssets);
% Estimated risk & return using EW port
% [EW_rsk, EW_ret] = estimatePortMoments(p, p.InitPort);

%% Set Efficient Frontier
% Set up portfolio constraints for fully-invested long-only portfolios
% (nonnegative weights that sum to 1)
p = setDefaultConstraints(p);
% Construct EFrontier
pwgt = estimateFrontier(p, 20);
% Vector of Risk & return across frontier
[frontier_rsk_vec, frontier_return_vec] = estimatePortMoments(p, pwgt);
%% Optimize for Sharpe Ratio
% Estimate portfolio wts with maximum SR objective
SR_wts = estimateMaxSharpeRatio(p,'Method','iterative','TolX', 1e-15);
SR_opt_wts = SR_wts';
% Estimate risk and return for optimal port
[MV_SR_sdv, MV_SR_return] = estimatePortMoments(p,SR_wts);

end





