function [simulationCellArray] = test_simulate_backtest(configurationMatrix, base_backtestObject)
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