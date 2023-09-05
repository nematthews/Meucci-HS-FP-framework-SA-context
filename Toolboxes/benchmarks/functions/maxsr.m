function [SR_opt_wts,MV_SR_sdv, MV_SR_return,frontier_rsk_vec, frontier_return_vec] = maxsr(AssetList,AssetMean,AssetCovar, Rfr,Wts_lb,Wts_up)
% Using MV asset allocation

% Mean-Variance Optimisation with a Maximum Sharpe Ratio objective function
% for portfolio weights allocation. Portfolio constraints use are that
% the portfolio is fully-invested long-only portfolios i.e.
% (nonnegative weights that sum to 1).

% INPUTS:
% AssetReturns: object consisting of J assets and their returns 
% over time period T
% (Type: timetable [T x J])

%% Set up data for MV SR portfolio object


%% Set up portfolio object
p = Portfolio('AssetList', AssetList, 'RiskFreeRate', Rfr);
p = setAssetMoments(p, AssetMean, AssetCovar);
p = setDefaultConstraints(p);

%% Set Weight Bounds (5 assets)
if nargin > 4
    p = setBounds(p, Wts_lb, Wts_up);
end
%% Set Efficient Frontier
% Set up portfolio constraints for fully-invested long-only portfolios
% (nonnegative weights that sum to 1)
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





