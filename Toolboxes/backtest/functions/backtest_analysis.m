function [hsfp_mu , hsfp_cov] = backtest_analysis(return_data,Window,Rfr,p_type)

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

%%

tickersToExtract = {'JALSHTR_Index', 'JIBA3M_Index','ALBTR_Index'};
BF_BH_TT = return_data(:, tickersToExtract);

[m,n]=size(return_data); % size of full data set
AssetList = return_data.Properties.VariableNames;
Rfr = 0.0056;
%###### initialize storage and inputs:

% ### 1. Equally Weighted (EW)###
Overlap_tsEW_Wts0 = eqweight(return_data); % [equity, cash, bonds] - balanced fund 60:40 
Overlap_tsEW_Wts = zeros(m,n); 
Overlap_tsEW_PRet = zeros(m,1); % storage for Portfolio Ret
% ### 2. SR Maximizing ###
Overlap_tsSR_Wts = zeros(m,n); % storage for Portfolio Weights
Overlap_tsSR_PRet = zeros(m,1); 
% ### 3. Balanced Fund Buy-Hold ###
Overlap_tsBH_Wts = [0.60,0.05,0.35]; % [equity, cash, bonds] - balanced fund 60:40 
% Overlap-window weights beginning of months
Overlap_tsBH_Wts0 = zeros(m,3);   % for t-th month
% Overlap-window weights end of months
Overlap_tsBH_WtsEnd = zeros(m,3); 
Overlap_tsBH_PRet = zeros(m,1);
% ### 4. HRP ###
Overlap_tsHRP_Wts = zeros(m,n);
Overlap_tsHRP_PRet = zeros(m,1); 
% ### 5. Balanced Fund Constant Mix (CM) ###
Overlap_tsCM_Wts0 = [0.60,0.05,0.35]; % [equity, cash, bonds] - balanced fund 60:40 
Overlap_tsCM_Wts = zeros(m,3); 
Overlap_tsCM_PRet = zeros(m,1); 


for t=Window:m-1
   
%%%%% need new stats of new window each time %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  m_t = mean(return_data{1+t-Window:t-1, :});
  cov_t = cov(return_data{1+t-Window:t-1,:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%% WEIGHTINGS %%%%%%
  % 1. EW
    % initialise wts as equally weighted
    Overlap_tsEW_Wts(t,:) = Overlap_tsEW_Wts0; % Constant Mix (CM)

  % 2. SR
    % initialise wts as equally weighted
    Overlap_tsSR_Wts(t,:) = maxsr(AssetList, m_t,cov_t, Rfr);

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
    Overlap_tsEW_PRet(t,:) = Overlap_tsEW_Wts(t,:) * transpose(return_data{t,:});
  % 2. SR
    Overlap_tsSR_PRet(t,:) =  Overlap_tsSR_Wts(t,:) * transpose(return_data{t,:});

  % 3. Balanced BH 
    Overlap_tsBH_PRet(t,:) = Overlap_tsBH_Wts0(t,:) * transpose(BF_BH_TT{t,:});
    
  % 4. HRP 
    Overlap_tsHRP_PRet(t,:) = Overlap_tsHRP_Wts(t,:) * transpose(return_data{t,:});

  % 5. Balanced CM 
    Overlap_tsCM_PRet(t,:) = Overlap_tsCM_Wts(t,:) * transpose(BF_BH_TT{t,:});

 %%%%% UPDATE BH %%%%%%%
 % Calc month end weight
 Overlap_tsBH_WtsEnd(t,:) = (BF_BH_TT{t,:}.*Overlap_tsBH_Wts0(t,:))+ Overlap_tsBH_Wts0(t,:);
 % Calc month (i+1) weights
 Overlap_tsBH_Wts0(t+1,:) = Overlap_tsBH_WtsEnd(t,:)/sum(Overlap_tsBH_WtsEnd(t,:));
end

end

