function [Opt_tsWts_cellarray,t] = backtest_analysis(backtest_object,Window)

% NOTE: This function is very use specific to this project. It was created
% to streamline the project code therefore does not generalise well.

%% INPUTS:
% returns_data - historical returns data for J assets
% (Type: TimeTable [T x J])

% window - window length using the same sample frequency as returns data
% (Type: scalar)

% Rfr - Risk free rate used for sharpe ratio optimisation usually the
% mean() of a cash asset calculated at the same frequency as returns.
% (Type: scalar)

% p_type - Method to apply with calculating HSFP flexible probabilities.
% (Type: char|str )

%% Sety up storage and initialise weights
returns_data = backtest_object.returns;

tickersToExtract = {'JALSHTR_Index', 'Cash','ALB_Index'};
BF_BH_TT = returns_data(:, tickersToExtract);

%###### initialize storage and inputs:
[m,n]=size(returns_data); % size of full data set
% Window = 36;
AssetList = returns_data.Properties.VariableNames;
% Rfr = mean(M_Ret_TT.JIBA3M_Index);
%###### initialize storage and inputs:

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

%% Use backtest objects for backtest moment calculations


%% Begin bakctest window shifts
for t=Window:m-1
    backtest_object.returns = returns_data{1+t-Window:t-1, :};

    % Need to calc pr for each window (but type stays the same)
    %  [hsfp_mu, hsfp_cov,p] = backtest_moments(backtest_object)
    %%%%% need new stats of new window each time %%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_t = exp(mean(log(returns_data{1+t-Window:t-1, :}+1)))-1;

    %   m_t = mean(M_Ret_TT{1+t-Window:t-1, :});
    cov_t = cov(returns_data{1+t-Window:t-1,:});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% WEIGHTINGS %%%%%%
    % 1. EW
    % initialise wts as equally weighted
    Overlap_tsEW_Wts(t,:) = Overlap_tsEW_Wts0; % Constant Mix (CM)

    % 2. SR
    % initialise wts as equally weighted
    %     Overlap_tsSR_Wts(t,:) = maxsr(AssetList, m_t,cov_t, M_Ret_TT.Cash(t));
    Overlap_tsSR_Wts(t,:) = maxsr(AssetList, m_t,cov_t, returns_data.Cash(t-1,:));

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

%% Trim resulting Weight window & store in cell array
% NOTE: each set of optimal weights have x number of columns based on numb
% assets. eg SR has 6 but BF have 3.

% Initialize the cell array
Opt_tsWts_cellarray = cell(5, 1);

% Assign the vectors to each cell in the cell array
Opt_tsWts_cellarray{1} = Overlap_tsEW_Wts(Window:t,:);
Opt_tsWts_cellarray{2} = Overlap_tsSR_Wts(Window:t,:);
Opt_tsWts_cellarray{3} = Overlap_tsBH_Wts0(Window:t,:);
Opt_tsWts_cellarray{4} = Overlap_tsHRP_Wts(Window:t,:);
Opt_tsWts_cellarray{5} = Overlap_tsCM_Wts(Window:t,:);

end

