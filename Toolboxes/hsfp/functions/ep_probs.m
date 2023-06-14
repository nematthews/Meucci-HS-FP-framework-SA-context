function [post_pr] = ep_probs(signal_series, alpha, z_target, prior)
% Flexible probabilities conditioned via entropy pooling

%%%%%%% ######## FIX INPUT TEXT TO OWN TEXT ###########

% INPUT
% Conditioner: [struct] with fields
% Series: [vector] (1 x t_) time series of the conditioner
% TargetValue: [vector] (1 x k_) target values for the conditioner
% Leeway: [scalar] (alpha) probability contained in the range, which is symmetric around the target value.
% prior: [vector] (1 x t_) prior set of probabilities

%% Crisp Prs
[p_cr] = cr_probs(signal_series, alpha, z_target);

% Crisp Mean & Sigma^2
cr_mu = sum(p_cr*signal_series);
cr_variance = sum(p_cr*(signal_series.^2))- cr_mu^2;

%% Entropy Pooling Setup %%

% Inequality constraints based on 2nd moments of crisp:
% z^2 <= (mu)^2 + varaince
a = (signal_series'.^2);   
b = (cr_mu^2)+cr_variance;   

% Equality constraints      
aeq = [signal_series';ones(1,length(signal_series))]; % Pr sum to 1
beq = [cr_mu;1];   % Constrain the 1st moments

post_pr = min_relative_entropy(prior', a,b,aeq,beq);
end