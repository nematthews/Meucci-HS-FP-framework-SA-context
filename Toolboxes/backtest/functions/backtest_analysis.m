function [Rolling_portfolioSet,Realised_tsPIndx,Realised_tsPRet,Opt_tsWts,cov_con_n,t,hsfp_Prs] = backtest_analysis(backtest_object,Window,reg_lambda)

% NOTE: This function is very use specific to this project. It was created
% to streamline the project code therefore does not generalise well.
% Function can also be made more streamline but currently coded in a
% practical quick manner.
%
%% Primary INPUT:
%% backtest_object - contains all data needed to calculate mu & sigma given
% a selected method.
% (Type: 1x1 struct with 6 fields)
%
% backtest_object = struct('returns',[],'Wts_lb',[],'Wts_ub',[],'signals',[],'method',[], 'parameters',[]);
%
% 1. backtest_object.returns - historical returns data for J assets
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
% 4. backtest_object.signals - time series of Q state signals.
% (Type: Timetable, [T x Q])
%
% 5. backtest_object.method - Method to calculate HSFP flexible probabilities.
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
% 6. backtest_object.parameters - contains fields for any parameters needed
% to calculate FProbs given an method that is selected.
% (Type: 1x1 struct with 8 fields)
% hsfp_parameters = struct('window',[],'tau',[],'alpha',[],'z_target',[],
% 'gamma',[],'h',[], 'tau_prior',[],'ensemble_wt_method',[]);
%
%% Additional INPUTS 
% Window - window length using the same sample frequency as returns data
% (Type: scalar)
%
% reg_lambda - regularisation parameter used to decrease noise within the
% covaraince matrix by means of adding a penalty term to the matrx. 
% (Type: double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example: backtest_object expected  (this case window needed for rolling_w)
% hsfp_parameters = struct('window',[],'tau',[],'alpha',[],'z_target',[],
% 'gamma',[],'h',[], 'tau_prior',[],'ensemble_wt_method',[]);
% hsfp_parameters.window = Window_len (NOTE: cannot exceed Window length
% specified in backtest_analysis.m Fn.)
% backtest_object = struct('returns',[],'Wts_lb',[],'Wts_ub',[],'signals',
% [],'method',[], 'parameters',[]);
% backtest_object.returns = returns_TT;
% backtest_object.Wts_lb = [ 0 0 0 0 0 ]; 
% backtest_object.Wts_ub = [ 1 1 1 1 0.05];
% backtest_object.signals = SIG_SMOOTHED_TT(:,'lagged_SACPIYOY_Index'); % NB: needs to be a timetable object
% backtest_object.method = 'rolling_w';
% backtest_object.parameters = hsfp_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NB NOTE: if backtest_object.method = 'rolling_w' ensure w_len < Window.

% Author: Nina Matthews (2023)

% $Revision: 1.1 $ $Date: 2023/05/09 19:09:01 $ $Author: Nina Matthews $


%% 1. Set up storage and initialise weights for Backtest %%%%%%%%%%%
portfolio_names = {'EW (CM)', 'MVSR max', 'BalFund (BH)', 'HRP', ...
    'BalFund (CM)', 'ALSI', 'ALBI', 'Cash'};

returns_data = backtest_object.returns;
signals_data = backtest_object.signals;

tickersToExtract = {'JALSHTR_Index', 'Cash','ALB_Index'};
BF_BH_TT = returns_data(:, tickersToExtract);
returns_data = removevars(returns_data,'JALSHTR_Index');

%###### initialize inputs:
[m,n]=size(returns_data); % size of full data set
% Window = 36;
AssetList = returns_data.Properties.VariableNames;
% Rfr = mean(returns_data.JIBA3M_Index);

%###### initialize storage:
% ### 0. covariance condition numbers
cov_con_n = zeros(m,1);
% ### 1. Equally Weighted (EW)###
Overlap_tsEW_Wts0 = eqweight(returns_data);
Overlap_tsEW_Wts = zeros(m,n);
Overlap_tsEW_PRet = zeros(m,1); % storage for Portfolio Ret
% ### 2. SR Maximizing ###
Overlap_tsSR_Wts = zeros(m,n); % storage for Portfolio Weights
Overlap_tsSR_PRet = zeros(m,1);
% ### 3. Balanced Fund Buy-Hold ###
Overlap_tsBH_Wts = [0.60,0.05,0.35]; % [equity, cash, bonds] - balanced fund 60:5:35
% Overlap-window weights beginning of months
Overlap_tsBH_Wts0 = zeros(m,3);   % for t-th month
% Overlap-window weights end of months
Overlap_tsBH_WtsEnd = zeros(m,3);
Overlap_tsBH_PRet = zeros(m,1);
% ### 4. HRP ###
Overlap_tsHRP_Wts = zeros(m,n);
Overlap_tsHRP_PRet = zeros(m,1);
% ### 5. Balanced Fund Constant Mix (CM) ###
Overlap_tsCM_Wts0 = [0.60,0.05,0.35]; % [equity, cash, bonds] - balanced fund 60:5:35
Overlap_tsCM_Wts = zeros(m,3);
Overlap_tsCM_PRet = zeros(m,1);

%% Set weight constraints for asset classes
if isempty(backtest_object.Wts_lb) 
    backtest_object.Wts_lb = zeros(1,n);
end

if isempty(backtest_object.Wts_ub) 
    backtest_object.Wts_ub = ones(1,n);
end

%% 2. Begin backtest window shifts %%%%%%%%%%%
for t=Window:m-1
    %%%%% need new stats of new window each time %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Instead of adding the conventional moments calc into
    % 'backtest_moments.m' we leave it independent so that the backtest still
    % functions if the file dependency fails for some reason.

    % Sets 'none' as default if backtest_object.method is not defined
    if strcmp(backtest_object.method, 'none') || isempty(backtest_object.method)
        % Geometric Mean:
        m_t = exp(mean(log(returns_data{1+t-Window:t-1, :}+1)))-1;
        % Arithmetic Covariance:
        cov_t = cov(returns_data{1+t-Window:t-1,:});
    else
        % Shift windows of returns and signals for HSFP each loop
        backtest_object.returns = returns_data(1+t-Window:t-1, :);
        backtest_object.signals = signals_data(1+t-Window:t-1, :);
        % Need to calc pr for each window (but type stays the same)
        [m_t, cov_t,hsfp_Prs] = backtest_moments(backtest_object);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Regularization parameter (lambda)
    % If reg_lambda is specified regularise cov matrix else skip:
    if nargin > 2 
    % Identity matrix of the same size as Cov
    n = size(cov_t,1);  % Assuming S is a square matrix
    I = eye(n);
    % Compute the regularized covariance matrix (R)
    cov_t = cov_t + reg_lambda * I;
    end

    % Checking condition number of cov at each t:
    cov_con_n(t,:) = cond(cov_t);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% WEIGHTINGS %%%%%%
    % 1. EW
    % initialise wts as equally weighted
    Overlap_tsEW_Wts(t,:) = Overlap_tsEW_Wts0; % Constant Mix (CM)

    % 2. SR
    % initialise wts as equally weighted
    %     Overlap_tsSR_Wts(t,:) = maxsr(AssetList, m_t,cov_t, returns_data.Cash(t));
    Overlap_tsSR_Wts(t,:) = maxsr(AssetList, m_t,cov_t, returns_data.Cash(t-1,:),backtest_object.Wts_lb,backtest_object.Wts_ub);

    % 3. Balanced BH
    % initialise wts as equally weighted
    Overlap_tsBH_Wts0(Window,:) = Overlap_tsBH_Wts; % initial weights

    % 4. HRP
    Overlap_tsHRP_Wts(t,:) = hrpestimate(cov_t)';

    % 5. Balanced CM
    % 3. Rebalance weighs
    Overlap_tsCM_Wts(t,:) = Overlap_tsCM_Wts0; % Constant Mix (CM)

    %%%%% PORT RETURNS %%%%%%
    % 1. EW
    Overlap_tsEW_PRet(t+1,:) = Overlap_tsEW_Wts(t,:) * transpose(returns_data{t+1,:});
    % 2. SR
    Overlap_tsSR_PRet(t+1,:) =  Overlap_tsSR_Wts(t,:) * transpose(returns_data{t+1,:});

    % 3. Balanced BH
    Overlap_tsBH_PRet(t+1,:) = Overlap_tsBH_Wts0(t,:) * transpose(BF_BH_TT{t+1,:});

    % 4. HRP
    Overlap_tsHRP_PRet(t+1,:) = Overlap_tsHRP_Wts(t,:) * transpose(returns_data{t+1,:});

    % 5. Balanced CM
    Overlap_tsCM_PRet(t+1,:) = Overlap_tsCM_Wts(t,:) * transpose(BF_BH_TT{t+1,:});

    %%%%% UPDATE BH %%%%%%%
    % Calc month end weight
    Overlap_tsBH_WtsEnd(t,:) = (BF_BH_TT{t,:}.*Overlap_tsBH_Wts0(t,:))+ Overlap_tsBH_Wts0(t,:);
    % Calc month (i+1) weights
    Overlap_tsBH_Wts0(t+1,:) = Overlap_tsBH_WtsEnd(t,:)/sum(Overlap_tsBH_WtsEnd(t,:));
end


%% 3. Alternative Benchmarks: ALSI & ALBI & Cash (JIBA3M) %%%%%%%%%%%
Overlap_tsALSI_PRet = table2array(BF_BH_TT(Window:t, 'JALSHTR_Index'));
% 0 to start of series to match length of return series from backtest loop
Overlap_tsALSI_PRet(1, :) = 0;
Overlap_tsALBI_PRet = table2array(returns_data(Window:t, 'ALB_Index'));
Overlap_tsALBI_PRet(1, :) = 0;
Overlap_tsCash_PRet = table2array(returns_data(Window:t, 'Cash'));
Overlap_tsCash_PRet(1, :) = 0;


%% 4. TRIM CONDITION NUMBERS:
% Trim resulting to window
cov_con_n = cov_con_n(Window:t,:);

%% 5. Package WEIGHTS for OUTPUT:
% Trim resulting to window & store in cell array %%%%%%%%%%%
% NOTE: each set of optimal weights have x number of columns based on numb
% assets. eg SR has 6 but BF have 3.

% Initialize the cell array
Opt_tsWts = cell(2, 5);

% Assign the vectors to each cell in the cell array
% Have portfolio Type as seperate arrays to identify cells
Opt_tsWts(1,:) = portfolio_names(1:5);
Opt_tsWts{2,1} = Overlap_tsEW_Wts(Window:t,:);
Opt_tsWts{2,2} = Overlap_tsSR_Wts(Window:t,:);
Opt_tsWts{2,3} = Overlap_tsBH_Wts0(Window:t,:);
Opt_tsWts{2,4} = Overlap_tsHRP_Wts(Window:t,:);
Opt_tsWts{2,5} = Overlap_tsCM_Wts(Window:t,:);

%% 6. Package RETURNS for OUTPUT:
% Trim resulting to window & store in cell array %%%%%%%%%%%
% Initialize the cell array
Realised_tsPRet = cell(2, 8);


% Assign the names to the first row of the cell array
Realised_tsPRet(1,:) = portfolio_names;
Realised_tsPRet{2,1} = Overlap_tsEW_PRet(Window:t,:);
Realised_tsPRet{2,2} = Overlap_tsSR_PRet(Window:t,:);
Realised_tsPRet{2,3} = Overlap_tsBH_PRet(Window:t,:);
Realised_tsPRet{2,4} = Overlap_tsHRP_PRet(Window:t,:);
Realised_tsPRet{2,5} = Overlap_tsCM_PRet(Window:t,:);
Realised_tsPRet{2,6} = Overlap_tsALSI_PRet;
Realised_tsPRet{2,7} = Overlap_tsALBI_PRet;
Realised_tsPRet{2,8} = Overlap_tsCash_PRet;


%% 7. Backtest Portfolio SHARPE RATIOS %%%%%%%%%%%
% ### 1. Equally Weighted (EW)###
SR_Overall_EW = sqrt(12)*((mean(Overlap_tsEW_PRet(Window:t,:))- ...
    mean(returns_data.Cash(Window:t,:)))/std(Overlap_tsEW_PRet(Window:t,:)));
% ### 2. SR Maximizing ###
SR_Overall_MVSR = sqrt(12)*((mean(Overlap_tsSR_PRet(Window:t,:))- ...
    mean(returns_data.Cash(Window:t,:)))/std(Overlap_tsSR_PRet(Window:t,:)));
% ### 3. Balanced Fund Buy-Hold ###
SR_Overall_BH = sqrt(12)*((mean(Overlap_tsBH_PRet(Window:t,:))- ...
    mean(returns_data.Cash(Window:t,:)))/std(Overlap_tsBH_PRet(Window:t,:)));
% ### 4. HRP ###
SR_Overall_HRP = sqrt(12)*((mean(Overlap_tsHRP_PRet(Window:t,:))- ...
    mean(returns_data.Cash(Window:t,:)))/std(Overlap_tsHRP_PRet(Window:t,:)));
% ### 5. Balanced Fund Constant Mix (CM) ###
SR_Overall_CM = sqrt(12)*((mean(Overlap_tsCM_PRet(Window:t,:))- ...
    mean(returns_data.Cash(Window:t,:)))/std(Overlap_tsCM_PRet(Window:t,:)));
% ### 6. ALSI - Equity Proxy ###
SR_Overall_ALSI = sqrt(12)*((mean(Overlap_tsALSI_PRet)- ...
    mean(returns_data.Cash(Window:t,:)))/std(Overlap_tsALSI_PRet));
% ### 7. ALBI - Bonds Proxy ###
SR_Overall_ALBI = sqrt(12)*((mean(Overlap_tsALBI_PRet)- ...
    mean(returns_data.Cash(Window:t,:)))/std(Overlap_tsALBI_PRet));
% ### 8. JIBA3M - Cash Proxy ###
SR_Overall_Cash = sqrt(12)*((mean(Overlap_tsCash_PRet)- ...
    mean(returns_data.Cash(Window:t,:)))/std(Overlap_tsCash_PRet));

%% 8. Portfolio PERFORMANCE: %%%%%%%%%%%
%  1. Calculate Geometrically Compounded Returns
%  2. Trim and store in cell array
% Preallocate for 8 variables (EW, SR, BH, HRP, CM, ALSI, ALBI, Cash)
Realised_tsPIndx = cell(2,8);
Realised_tsPIndx(1,:) = portfolio_names;

% Define the list of port types and their corresponding names
variables = {
    'EW', 'SR', 'BH', 'HRP', 'CM', 'ALSI', 'ALBI', 'Cash'
    };

% Loop through the ports, calculate the indices, and store
for i = 1:length(variables)
    % Calculate the decimal return and cumulative prod
    decimal_ret = eval(['Overlap_ts' variables{i} '_PRet(1:end,:) + 1']);
    PIndx = cumprod(decimal_ret);

    % Store the results in the cell array
    Realised_tsPIndx{2,i} = PIndx;
end


%% 9. Create Timetable object to store cummulative returns in %%%%%%%%%%%
Rolling_portfolioSet = timetable(returns_data.Time(Window:t,:), ...
    Realised_tsPIndx{2,1}(Window:t,:), ...
    Realised_tsPIndx{2,2}(Window:t,:), ...
    Realised_tsPIndx{2,3}(Window:t,:), ...
    Realised_tsPIndx{2,4}(Window:t,:), ...
    Realised_tsPIndx{2,5}(Window:t,:), ...
    'VariableNames',{'EW (CM)', ...
    'MVSR max', ...
    'BalFund (BH)', ...
    'HRP', ...
    'BalFund (CM)'});
%%% Adding in the additional benchmark series
Rolling_portfolioSet=[Rolling_portfolioSet timetable(returns_data.Time(Window:t,:), ...
    Realised_tsPIndx{2,6},Realised_tsPIndx{2,7}, ...
    Realised_tsPIndx{2,8} , ...
    'VariableNames',{'ALSI', ...
    'ALBI', ...
    'Cash'})];


%% 10. Create Timetable with SRs to store cummulative returns in %%%%%%%%%%%
Rolling_portfolioSet_SR = timetable(returns_data.Time(Window:t,:), ...
    Realised_tsPIndx{2,1}(Window:t,:), ...
    Realised_tsPIndx{2,2}(Window:t,:), ...
    Realised_tsPIndx{2,3}(Window:t,:), ...
    Realised_tsPIndx{2,4}(Window:t,:), ...
    Realised_tsPIndx{2,5}(Window:t,:), ...
    'VariableNames',{['EW (CM)  (SR = ' num2str(round(SR_Overall_EW,2)) ')'], ...
    ['MVSR max  (SR = ' num2str(round(SR_Overall_MVSR,2)) ')'], ...
    ['BalFund (BH)   (SR = ' num2str(round(SR_Overall_BH,2)) ')'], ...
    ['HRP   (SR = ' num2str(round(SR_Overall_HRP,2)) ')'], ...
    ['BalFund (CM)    (SR = ' num2str(round(SR_Overall_CM,2)) ')']});
%%% Adding in the additional benchmark series
Rolling_portfolioSet_SR=[Rolling_portfolioSet_SR timetable(returns_data.Time(Window:t,:), ...
    Realised_tsPIndx{2,6},Realised_tsPIndx{2,7}, ...
    Realised_tsPIndx{2,8} , ...
    'VariableNames',{['ALSI  (SR = ' num2str(round(SR_Overall_ALSI,2)) ')'], ...
    ['ALBI  (SR = ' num2str(round(SR_Overall_ALBI,2)) ')'], ...
    ['Cash (JIBA3M)  (SR = ' num2str(round(SR_Overall_Cash,2)) ')']})];
end

