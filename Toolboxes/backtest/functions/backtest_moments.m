function [hsfp_mu, hsfp_cov,p] = backtest_moments(backtest_object)

% NOTE: This function is very use specific to this project. It was created
% to streamline the project code therefore does not generalise well.
%
%% INPUTS:
%
%% backtest_object - contains all data needed to calculate mu & sigma given
% a selected method.
% (Type: 1x1 struct with 6 fields)
%
% backtest_object = struct('Returns',[],'Wts_lb',[],'Wts_ub',[],'Signals',[],'Methods',[], 'HSFPparameters',[]);
%
% 1. backtest_object.Returns - historical Returns data for J assets
% (Type: TimeTable [T x J])
%
% 2. backtest_object.Wts_lb - vector of lower bounds for the asset weights 
% in the portfolios. (Auxiliary (used elsewhere in backtest_analysis.m))
% (Type: double [1 x J])
%
% 3. backtest_object.Wts_ub - vector of upper bounds for the asset weights
% in the portfolios. (Auxiliary (used elsewhere in backtest_analysis.m))
% (Type: double [1 x J])
%
% 4. backtest_object.Signals - time series of Q state Signals.
% (Type: Timetable, [T x Q])
%
% 5. backtest_object.Methods - Method to calculate HSFP flexible probabilities.
% (Type: char|str )
% Options:   - none (default) NOTE: this is not defined below as it is
%                                   only used in backtest fn to rather use 
%                                   normal mean and cov. 
%            - rolling_w 
%            - exp_decay
%            - crisp
%            - kernel
%            - e_pooling
%            - ew_ensemble
%            - cb_ensemble
%
% 6. backtest_object.HSFPparameters - contains fields for any parameters needed
% to calculate FProbs given an Methods that is selected.
% (Type: 1x1 struct with 8 fields)
% hsfp_parameters = struct('RollWindow',[],'Tau',[],'Alpha',[],'Z_target',[],
% 'Gamma',[],'h',[], 'Tau_prior',[],'Ensemble_wt_method',[]);
%
%% Auxiliary INPUTS (used elsewhere in backtest_analysis.m):
% Window - RollWindow length using the same sample frequency as Returns data
% (Type: scalar)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example: backtest_object expected  (this case window needed for rolling_w)
% hsfp_parameters = struct('RollWindow',[],'Tau',[],'Alpha',[],'Z_target',[],
% 'Gamma',[],'h',[], 'Tau_prior',[],'Ensemble_wt_method',[]);
% hsfp_parameters.RollWindow = Window_len (NOTE: cannot exceed Window length
% specified in backtest_analysis.m Fn.)
% backtest_object = struct('Returns',[],'Wts_lb',[],'Wts_ub',[],'Signals',
% [],'method',[], 'HSFPparameters',[]);
% backtest_object.Returns = Returns_TT;
% backtest_object.Wts_lb = [ 0 0 0 0 0 ]; 
% backtest_object.Wts_ub = [ 1 1 1 1 0.05];
% backtest_object.Signals = SIG_SMOOTHED_TT(:,'lagged_SACPIYOY_Index'); % NB: needs to be a timetable object
% backtest_object.method = 'rolling_w';
% backtest_object.HSFPparameters = hsfp_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required parameters per method:
%            - rolling_w    (RollWindow)
%            - exp_decay    (Tau)
%            - crisp        (Alpha, Z_target)
%            - kernel       (h, Gamma, Z_target)
%            - e_pooling    (Alpha, Tau_prior, Z_target)
%            - ew_ensemble  (Alpha, Tau_prior, Z_target)
%            - cb_ensemble  (Alpha, Tau_prior, Z_target,Ensemble_wt_method)


% Author: Nina Matthews (2023)

% $Revision: 1.5 $ $Date: 2023/05/09 19:09:01 $ $Author: Nina Matthews $

%% Selection of HS-FP method for calculating Flexible Probabilites
if strcmp(backtest_object.Method, 'rolling_w')
    % Rolling Window Prs
    p = rw_probs(backtest_object.Signals,backtest_object.HSFPparameters.RollWindow);

elseif strcmp(backtest_object.Method, 'exp_decay')
    p = es_probs(backtest_object.Signals, backtest_object.HSFPparameters.Tau);

elseif strcmp(backtest_object.Method, 'crisp')
    p = cr_probs(backtest_object.Signals, backtest_object.HSFPparameters.Alpha, ...
        backtest_object.HSFPparameters.Z_target);

elseif strcmp(backtest_object.Method, 'kernel')
    p = kl_probs(backtest_object.Signals, backtest_object.HSFPparameters.h, ...
        backtest_object.HSFPparameters.Gamma,backtest_object.HSFPparameters.Z_target);

elseif strcmp(backtest_object.Method, 'e_pooling')
    p = ep_probs(backtest_object.Signals, backtest_object.HSFPparameters.Alpha, ...
        backtest_object.HSFPparameters.Tau_prior,backtest_object.HSFPparameters.Z_target);

elseif strcmp(backtest_object.Method, 'ew_ensemble')
    p = ew_ensemble(backtest_object.Signals, backtest_object.HSFPparameters.Alpha, ...
        backtest_object.HSFPparameters.Tau_prior,backtest_object.HSFPparameters.Z_target);

elseif strcmp(backtest_object.Method, 'cb_ensemble')
    p = cb_ensemble(backtest_object.Signals, backtest_object.HSFPparameters.Alpha, ...
        backtest_object.HSFPparameters.Tau_prior,backtest_object.HSFPparameters.Z_target, ...
        backtest_object.HSFPparameters.Ensemble_wt_method);
end

[hsfp_mu , hsfp_cov] = hsfp_moments(p,backtest_object.Returns);

% Calc mean and var, add condition for using conventional calc
end



