function [data_simulationCellArray,parameter_simulationCellArray] = para_simulate_backtest(data_configurationMatrix,parameter_configurationMatrix,base_backtestObject)
% Generates multiple backtest simulations when given a configuration matrix
% that contains multiple combinations of varying hyper parameters. The
% configuartion matrix can be constructed manually based on a users choice
% of specific configuartions that they would like to assess or it can be
% constructed using a range of values for each hyper parameter and finding
% all j possible combinations of those values. The output of the function
% results in a cell array object with j different backtestedPortfolio
% objects.
%
%% INPUT:
%
% data_configurationMatrix - matrix of all possible (j) combinations
% of all (h) varying hyper-parameters needed for data processing testing.
% Should include combiniations of parameters such as cash contraints, winsoring
% standard deviations and regularisation lambdas.
% Where j = length of range of h^h.
%
% parameter_configurationMatrix - matrix of all possible (j) combinations
% of all (h) varying HSFP hyper-parameters needed for sensitivity anaylsis of the HSFP methods.
% Should include combiniations of parameters such as rolling window length,
% tau, gamma, alpha etc. (see HSFPparameter class definition for details).
% Where j = length of range of h^h.
%
% Eg: with 3 hyper parameters with a range of 15 values, j = 15^3.
% (type: array double, [j x h] | double)
%
% base_backtestObject - this backtest object has the desired constant
% parameters defined as the objects properties.
% (type: backtestedPortfolio class object)


% See also backtestedPortfolios.m

% Author: Nina Matthews (2023)


% $Revision: 1.2 $ $Date: 2023/10/26 10:23:41 $ $Author: Nina Matthews $

%% Run simulations
%% Apply data processing configurations
if ~isempty(data_configurationMatrix)
    %%% Create storage for simulation objects
    num_data_configs = size(data_configurationMatrix, 1);
    data_simulationCellArray = cell(1, num_data_configs);

    % Create a cell array to store the test objects
    data_test_class1Array = cell(1, num_data_configs);

    % Populate the cell array with initial objects
    parfor config = 1:num_data_configs
        iterativeClass = copy(base_backtestObject);  % Use the copy method
        iterativeClass.CashConstriant = data_configurationMatrix(config, 1);
        iterativeClass.WinsorStd = data_configurationMatrix(config, 2);
        iterativeClass.RegLambda = data_configurationMatrix(config, 3);
        data_test_class1Array{config} = iterativeClass;
    end

    % Use parfor to run backtests
    parfor config = 1:num_data_configs
        data_simulationCellArray{config} = OOPbacktest_analysis(data_test_class1Array{config});
    end
    % will add an if statement where both data and parameter simulations
    % take place
    parameter_simulationCellArray = [];
end

%% Apply parameter configurations
if  ~isempty(parameter_configurationMatrix)
    %%% Create storage for simulation objects
    parameter_num_configs = size(parameter_configurationMatrix, 1);
    parameter_simulationCellArray = cell(1, parameter_num_configs);

    % Create a cell array to store the test objects
    parameter_test_class1Array = cell(1, parameter_num_configs);

    % Populate the cell array with initial objects
    parfor config = 1:parameter_num_configs
        iterativeClass = copy(base_backtestObject);  % Use the copy method
        iterativeClass.HSFPparameters.RollWindow  = parameter_configurationMatrix(config, 1);
        % iterativeClass.HSFPparameters.Tau = parameter_configurationMatrix(config, 2);
        % iterativeClass.HSFPparameters.Z_target = parameter_configurationMatrix(config, 3);
        % iterativeClass.HSFPparameters.Alpha = parameter_configurationMatrix(config, 4);
        % iterativeClass.HSFPparameters.h = parameter_configurationMatrix(config, 5);
        % iterativeClass.HSFPparameters.Tau_prior = parameter_configurationMatrix(config, 6);
        % iterativeClass.HSFPparameters.Gamma = parameter_configurationMatrix(config, 7);
        parameter_test_class1Array{config} = iterativeClass;
    end

    % Use parfor to run backtests
    parfor config = 1:parameter_num_configs
        parameter_simulationCellArray{config} = OOPbacktest_analysis(parameter_test_class1Array{config});
    end

    data_simulationCellArray = [];
end
end




