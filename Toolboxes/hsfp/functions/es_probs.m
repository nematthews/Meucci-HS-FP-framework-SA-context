function p_es = es_probs(Sig_timeTable, tau)
% exponential smoothing
% es_probs facilitates the calculation of exponentially smoothed 
% time-conditioned probabilities used in the HS-FP framework 
% developed by Meucci (2010). 
% 
%% INPUTS:
% Sig_timeTable - timeseries of state signals 
% (Type: TimeTable object)
%
% Tau - half-life decay (time it would take for Pr to decay to half
% the eval of T.)
% (Type: double)


% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/02/20 16:10:46 $ $Author: Nina Matthews $

%% 
% Final point of Timetable 
T = height(Sig_timeTable);

% check default value: tau = 24 mnths (2 yrs)
if nargin < 2 || isempty(tau)
    tau = 24; 
end
   
% exp decay fn
p = exp(-(log(2)/tau)*(T-(1:T)'));
% Rescale
p_es = (p/sum(p))';

end