function [hsfp_mu, hsfp_cov,p] = backtest_moments(backtest_object,Window,t)

% NOTE: This function is very use specific to this project. It was created
% to streamline the project code therefore does not generalise well.

%% INPUTS:

% backtest_object - contains all data needed to calculate mu & sigma given
% a selected method.
% (Type: 1x1 struct with 4 fields)
% backtest_object = struct('returns',[],'signals',[],'method',[], 'parameters',[]);

% 1. backtest_object.returns - historical returns data for J assets
% (Type: TimeTable [T x J])

% 2. backtest_object.signals - time series of Q state signals.
% (Type: Timetable, [T x Q])

% 3. backtest_object.method - Method to calculate HSFP flexible probabilities.
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

% 4. backtest_object.parameters - contains fields for any parameters needed
% to calculate FProbs given an method that is selected.
% (Type: 1x1 struct with 9 fields)
% hsfp_parameters = struct('window',[],'tau',[],'alpha',[],'z_target',[],
% 'gamma',[],'h',[], 'EPool_prior',[], 'z_target_method', [],'ensemble_wt_method',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example: backtest_object expected  (this case window needed for rolling_w)
% hsfp_parameters = struct('window',[],'tau',[],'alpha',[],'z_target',[],
% 'gamma',[],'h',[], 'EPool_prior',[], 'z_target_method', [],'ensemble_wt_method',[]);
% hsfp_parameters.window = Window_len (NOTE: cannot exceed Window length
% specified in backtest_analysis.m Fn.)
% backtest_object = struct('returns',[],'signals',[],'method',[], 'parameters',[]);
% backtest_object.returns = returns_TT;
% backtest_object.signals = SIG_SMOOTHED_TT(:,'lagged_SACPIYOY_Index'); % NB: needs to be a timetable object
% backtest_object.method = 'rolling_w';
% backtest_object.parameters = hsfp_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Selection of HS-FP method for calculating Flexible Probabilites
if strcmp(backtest_object.method, 'rolling_w')
    % Rolling window Prs
    p = rw_probs(backtest_object.signals,backtest_object.parameters.window);

elseif strcmp(backtest_object.method, 'exp_decay')
    p = es_probs(backtest_object.signals, backtest_object.parameters.tau);

elseif strcmp(backtest_object.method, 'crisp')
    p = cr_probs(backtest_object.signals, backtest_object.parameters.alpha, ...
        backtest_object.parameters.z_target);

elseif strcmp(backtest_object.method, 'kernel')
    p = kl_probs(backtest_object.signals, backtest_object.parameters.h, ...
        backtest_object.parameters.gamma,backtest_object.parameters.z_target);

elseif strcmp(backtest_object.method, 'e_pooling')
    p = ep_probs(backtest_object.signals, backtest_object.parameters.alpha, ...
        backtest_object.parameters.z_target,backtest_object.parameters.EPool_prior);

elseif strcmp(backtest_object.method, 'ew_ensemble')
    p = ew_ensemble(backtest_object.signals, backtest_object.parameters.alpha, ...
        backtest_object.parameters.EPool_prior,backtest_object.parameters.z_target_method);

elseif strcmp(backtest_object.method, 'cb_ensemble')
    p = cb_ensemble(backtest_object.signals, backtest_object.parameters.alpha, ...
        backtest_object.parameters.EPool_prior,backtest_object.parameters.z_target_method);
end

[hsfp_mu , hsfp_cov] = hsfp_moments(p,backtest_object.returns);

% Calc mean and var, add condition for using conventional calc
end



