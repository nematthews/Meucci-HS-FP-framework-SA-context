function [x,retSR, rskSR] = maxsr(Returns, Rfr,lb,ub)

% Using MV asset allocation
%
% Mean-Variance Optimisation with a Maximum Sharpe Ratio objective function
% for portfolio weights allocation. Portfolio constraints use are that
% the portfolio is fully-invested long-only portfolios i.e.
% (nonnegative weights that sum to 1).
%
% NOTE: makes use of geometric means not arithmetic
%
%% INPUTS:
% AssetReturns: object consisting of J assets and their returns 
% over time period T
% (Type: timetable [T x J])

% Author: Nina Matthews (2023)

% $Revision: 1.1 $ $Date: 2023/05/09 19:09:01 $ $Author: Nina Matthews $

% risk free rate
RFR = Rfr;
% monthy sampled arithmetic mean of simple returns
mu = geo_ave(Returns);
% average risk free rate (when to use geometric average)
ERFR = geo_ave(RFR); 
% covariance
Sigma = cov(Returns);
% Equally weighted portfolio
x0 = ones(size(mu))/length(mu);
%initialise the weights
Wts = NaN(1,length(mu));
% objective function to maximise the SR
fn0 = @(x) (-(x*mu' - ERFR)/sqrt(x*Sigma*x'));
% % % No Short-selling (upper and lower bounds)
% % ub = ones(length(mu),1);
% % lb = zeros(length(mu),1);
% Equality constraint (fully invested)
Aeq = ones(1,length(mu));
beq = 1;
% Use SQP to solve for the tangency portfolio
options = optimoptions(@fmincon,'Algorithm','sqp','OptimalityTolerance',1e-8,'Display','off');
x = fmincon(fn0,x0,[],[],Aeq,beq,lb,ub,[],options); 
% plot the maximum SR portfolio on EF
retSR = x*mu';
rskSR = sqrt(x*Sigma*x');
% clf; portfolioexamples_plot('Efficient Frontier (QP)', ...
    % {'line', rsk, ret}, ...
    % {'scatter',rskSR,retSR,{'Maximum Sharpe ratio'},'r'}, ...
	% {'scatter', sqrt(diag(q.AssetCovar)), q.AssetMean, q.AssetList, '.r'});

end 