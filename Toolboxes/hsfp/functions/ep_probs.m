function [post_pr] = ep_probs(signal_series, alpha, z_target, prior)
% Flexible probabilities conditioned via entropy pooling

% Allows for Time & State conditioning of FPs by begining
% with an exponential decay prior and uses the Crisp Prs & its moments to 
% for condition. Due to the Entropy pooling we have a resulting mixture
% of the kernel approach and the exponential decay that can
% switch between Gaussian kernel and exponential kernel. 

% Calls on min_relative_entropy.m to calculate the posterior probabilites
% when the Crisp Moments are inputs as the desired constraints. 

% INPUT:
% singnal_series - smoothed & standardised timeseries for a single state var
% (type: array double, [T x 1])

% alpha - range of probability for bandwidth
% (type: double)

% z_target - target value specific to the signal_series
% (type: double)

% prior - desired prior probability distribution (usually exp smoothed Prs)
% (type: array double [T x 1]) 

%% Crisp Prs
p_cr = cr_probs(signal_series, alpha, z_target);
p_cr(p_cr==0)=10^-20;


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