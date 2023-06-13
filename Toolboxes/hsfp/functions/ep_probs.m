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

% Crisp Mean & Sigma
cr_mu = sum(p_cr*signal_series);
cr_sigma = sum(p_cr*(signal_series.^2))- cr_mu^2;

%% Entropy Pooling Setup %%
% Inequality constraints based on crisp
% 1st Moment to match: mu
a = (signal_series'.^2);   
% 2nd Moment to match: sigma
b = (cr_mu^2)+cr_sigma;  

% Equality constraints      
aeq = [signal_series';ones(1,length(signal_series))];
beq = [cr_mu;1];

post_pr = MinRelEntFP_001(prior', a,b,aeq,beq);
end