function [simulationCellArray] = simulate_backtest(configurationMatrix,base_backtestObject)
% Generates multiple backtest simulations when given a configuration matrix
% that contains multiple combinations of varying hyper parameters. The
% configuartion matrix can be constructed manually based on a users choice
% of specific configuartions that they would like to assess or it can be
% constructed using a range of values for each hyper parameter and finding
% all j possible combinations of those values. The out put of the function
% results in a cell array object with j different backtestedPortfolio
% objects. 
%
%% INPUT:
%
% configurationMatrix - matrix of all possible (j) combinations
% of all (h) varying hyper-parameters. Where j = length of range of h^h.
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
%%% Create storage for simulation objects
    num_configs = size(configurationMatrix, 1);
    simulationCellArray = cell(1, num_configs);

    % Create a cell array to store the test objects
    test_class1Array = cell(1, num_configs);

    % Populate the cell array with initial objects
    parfor config = 1:num_configs
        iterativeClass = copy(base_backtestObject);  % Use the copy method
        iterativeClass.CashConstriant = configurationMatrix(config, 1);
        iterativeClass.WinsorStd = configurationMatrix(config, 2);
        iterativeClass.RegLambda = configurationMatrix(config, 3);
        test_class1Array{config} = iterativeClass;
    end

    % Use parfor to run backtests
    parfor config = 1:num_configs
        simulationCellArray{config} = OOPbacktest_analysis(test_class1Array{config});
    end

end



