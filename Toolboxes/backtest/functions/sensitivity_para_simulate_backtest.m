function [parameter_sim_CellArray] = sensitivity_para_simulate_backtest(base_backtestObject, num_parameters, signal_series)
% Generates multiple backtest simulations specifically to simulate different
% HSFP methods when given a base backtestedPortfolio object, a number of
% values wanted to generate a range of values per hyper-parameter and a set
% of state signals to use in the HSFP framework.
% From the inputs, differrent configuration matrices are generated per
% method based on its specific hyper-parameters as well as the signals
% provided. Switch cases are based on specified method type in based
% backtestedPortfolio.
%
% NOTE: If data processing configurations want to be tested, first apply
% this using para_simulate_backtest.m function. Once desired data processes
% have been chosen based on DSR results, make sure to specify those chosen
% prcesses as the initial properties of the base backtestedPorflio object
% passed into this function.
%
% ** Step 1 **: generate parameter_configurationMatrix -
% matrix of all possible (j) combinations of all (h) varying HSFP
% hyper-parameters needed for sensitivity anaylsis of the HSFP methods.
% Should include combiniations of parameters such as rolling window length,
% tau, gamma, alpha etc. (see HSFPparameter class definition for details).
% Where j = length of range of h^h.
%
%% INPUT:
% base_backtestObject - this backtest object has the desired constant
% parameters defined as the objects properties.
% (type: backtestedPortfolio class object)
%
% signal_series - HSFP signals on which to condition. All signals provided will
% be used to generate alternative configuartions to be simulated for each
% method.
% (type: Timetable, [T x Q])
%
% num_parameters - number of different equally spaced values to generate
% per hyper-parameter. Eg: 3 hyper parameters with a range of 5 values, j = 5^3.
% (type: double)
%
%% OUTPUT:
% parameter_simulationCellArray - cell array object with M different
% backtestedPortfolio objects per cell. Each cell object is speficic to a
% unique configuration of signals and parameters.
% (type: cell array, {1 x number of configurations})
%
%% Options for each method: and required parameters
%            - none (default)
% Time & Single State Conditioning Methods
%            - rolling_w:   (window)
%            - exp_decay:    (tau)
%            - crisp:       (alpha, z_target [can be 'mean','latest' or scalar])
%            - kernel:      (h, gamma, z_target [can be 'mean','latest' or scalar])
%            - e_pooling:   (alpha, tau_prior, z_target [can be 'mean','latest' or scalar])
% Multi-State & Time Conditioning Methods
%            - ew_ensemble: (alpha, tau_prior, z_target [can ONLY be 'mean' or 'latest'])
%            - cb_ensemble: (alpha, tau_prior, z_target [can ONLY be 'mean' or 'latest'] ,ensemble_wt_method)
%
%
%% See also backtestedPortfolios.m, para_simulate_backtest.m

% Author: Nina Matthews (2023)


% $Revision: 1.2 $ $Date: 2023/12/01 10:23:41 $ $Author: Nina Matthews $

%% Setup simulation objects
% if isa(signal_series,'timetable')
%     signal_series = table2array(signal_series);
% end

num_signals = width(signal_series);
singal_index = 1:num_signals;

%% Method switch changes
switch base_backtestObject.Method

    case 'rolling_w'
        %% 1. Rolling Window Time conditioned %%%%%%%%%%%%%%%%%%%%%
        %%%% rolling_w:   (window) %%%%
        RollWindow_range = round(linspace(3, 35, num_parameters)); % Double. Rounded to allow for rolling window
        parameter_configMatrix = RollWindow_range';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%% Apply configuration in simulations %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Create storage for simulation objects
        parameter_num_configs = size(parameter_configMatrix, 1);
        parameter_sim_CellArray = cell(1, parameter_num_configs);

        % Create a cell array to store the test objects
        parameter_test_class1Array = cell(1, parameter_num_configs);

        % Populate the cell array with initial objects
        parfor config = 1:parameter_num_configs
            iterativeClass = copy(base_backtestObject);  % Use the copy method
            iterativeClass.HSFPparameters.RollWindow  = parameter_configMatrix(config, 1);
            % assign object
            parameter_test_class1Array{config} = iterativeClass;
        end

        % Use parfor to run backtests
        parfor config = 1:parameter_num_configs
            parameter_sim_CellArray{config} = OOPbacktest_analysis(parameter_test_class1Array{config});
        end


    case 'exp_decay'
        %% 2. Exponential Decay Time conditioned %%%%%%%%%%%%%%%%%%%%%
         %%%% exp_decay: (tau) %%%%
        Tau_range = linspace(3, 36, num_parameters); % Double
        parameter_configMatrix = Tau_range';

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%% Apply configuration in simulations %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% Create storage for simulation objects
        parameter_num_configs = size(parameter_configMatrix, 1);
        parameter_sim_CellArray = cell(1, parameter_num_configs);

        % Create a cell array to store the test objects
        parameter_test_class1Array = cell(1, parameter_num_configs);

        % Populate the cell array with initial objects
        parfor config = 1:parameter_num_configs
            iterativeClass = copy(base_backtestObject);  % Use the copy method
            iterativeClass.HSFPparameters.Tau  = parameter_configMatrix(config, 1);
            % assign object
            parameter_test_class1Array{config} = iterativeClass;
        end

        % Use parfor to run backtests
        parfor config = 1:parameter_num_configs
            parameter_sim_CellArray{config} = OOPbacktest_analysis(parameter_test_class1Array{config});
        end

    case 'crisp'
        %% 3. Crisp State conditioned %%%%%%%%%%%%%%%%%%%%%
        %%% Parameters need:
        %%%% crisp: (alpha, z_target [can be 'mean','latest' or scalar]) %%%%

        % storage array
        parameter_sim_CellArray = cell(1,1);

        Alpha_range = linspace(0.1, 1, num_parameters); % Double (needs to be at least 0.1 to allow convergence in MAXSR)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Loop through all signals %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for signal = 1:num_signals

            %%% Generate configurations %%%%

            % get specific Z_target per signal based on bounds of signal
            [Sig_min,Sig_max] = bounds(signal_series{:,signal});
            Z_target_range = linspace(Sig_min, Sig_max, num_parameters); % Generated based on signal range (range between min and max)

            hyperparameter_set = [Z_target_range; Alpha_range];

            % Initialize an empty matrix to store configurations
            num_configs = size(hyperparameter_set, 2)^size(hyperparameter_set, 1);
            parameter_configMatrix = zeros(num_configs, size(hyperparameter_set, 1));

            % Nested loops to generate all configurations based on the size of hyperparameter_set
            col = 1;
            for config1 = 1:size(hyperparameter_set, 2)
                for config2 = 1:size(hyperparameter_set, 2)
                    % Store the current combination of values in the configurations matrix
                    parameter_configMatrix(col, :) = [hyperparameter_set(1, config1); hyperparameter_set(2, config2)];
                    col = col + 1;
                end
            end

            base_backtestObject.Signals = signal_series(:,signal);

            %%%%%%%%%%%%%%% Apply configuration in simulations %%%%%%%

            %%% Create storage for simulation objects
            parameter_num_configs = size(parameter_configMatrix, 1);
            parameter_simulationCellArray = cell(1, parameter_num_configs);

            % Create a cell array to store the test objects
            parameter_test_class1Array = cell(1, parameter_num_configs);

            % Populate the cell array with initial objects
            parfor config = 1:parameter_num_configs
                iterativeClass = copy(base_backtestObject);  % Use the copy method
                iterativeClass.HSFPparameters.Z_target = parameter_configMatrix(config, 1);
                iterativeClass.HSFPparameters.Alpha = parameter_configMatrix(config, 2);
                % assign object
                parameter_test_class1Array{config} = iterativeClass;
            end

            % Use parfor to run backtests
            parfor config = 1:parameter_num_configs
                parameter_simulationCellArray{config} = OOPbacktest_analysis(parameter_test_class1Array{config});
            end

            % Accumulate results each signal loop
            parameter_sim_CellArray = horzcat(parameter_sim_CellArray,parameter_simulationCellArray);
        end

        % Remove the first cell
        parameter_sim_CellArray(1) = [];

    case 'kernel'
        %% 4. Kernel State conditioned %%%%%%%%%%%%%%%%%%%%%
        %%%% kernel: (h, gamma, z_target [can be 'mean','latest' or scalar]) %%%%

        % storage array
        parameter_sim_CellArray = cell(1,1);

        h_range = linspace(0.05, 1, num_parameters); % Double

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Loop through all signals %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for signal = 1:num_signals


            %%% Generate configurations %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % get specific Z_target per signal based on bounds of signal
            [Sig_min,Sig_max] = bounds(signal_series{:,signal});
            Z_target_range = linspace(Sig_min, Sig_max, num_parameters); % Generated based on signal range (range between min and max)

            hyperparameter_set = [h_range; Z_target_range];


            % Initialize an empty matrix to store configurations
            num_configs = size(hyperparameter_set, 2)^size(hyperparameter_set, 1);
            parameter_configMatrix = zeros(num_configs, size(hyperparameter_set, 1));

            % Nested loops to generate all configurations based on the size of hyperparameter_set
            col = 1;
            for config1 = 1:size(hyperparameter_set, 2)
                for config2 = 1:size(hyperparameter_set, 2)
                    % Store the current combination of values in the configurations matrix
                    parameter_configMatrix(col, :) = [hyperparameter_set(1, config1); hyperparameter_set(2, config2)];
                    col = col + 1;
                end
            end

            % Account for the binary option for gamma of kernel:
            % gamma = 1 is exponential
            % gamma = 2 is Gaussian
            configurations_gamma1 = [parameter_configMatrix, ones(size(parameter_configMatrix, 1), 1)];
            configurations_gamma2 = [parameter_configMatrix, 2*ones(size(parameter_configMatrix, 1), 1)];
            parameter_configMatrix = [configurations_gamma1; configurations_gamma2];

            % update signal per loop
            base_backtestObject.Signals = signal_series(:,signal);


            %%%%%%%%%%%%%%% Apply configuration in simulations %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% Create storage for simulation objects
            parameter_num_configs = size(parameter_configMatrix, 1);
            parameter_simulationCellArray = cell(1, parameter_num_configs);

            % Create a cell array to store the test objects
            parameter_test_class1Array = cell(1, parameter_num_configs);

            % Populate the cell array with initial objects
            parfor config = 1:parameter_num_configs
                iterativeClass = copy(base_backtestObject);  % Use the copy method
                iterativeClass.HSFPparameters.h = parameter_configMatrix(config, 1);
                iterativeClass.HSFPparameters.Z_target = parameter_configMatrix(config, 2);
                iterativeClass.HSFPparameters.Gamma = parameter_configMatrix(config, 3);
                % assign object
                parameter_test_class1Array{config} = iterativeClass;
            end

            % Use parfor to run backtests
            parfor config = 1:parameter_num_configs
                parameter_simulationCellArray{config} = OOPbacktest_analysis(parameter_test_class1Array{config});
            end
            % Accumulate results each signal loop
            parameter_sim_CellArray = horzcat(parameter_sim_CellArray,parameter_simulationCellArray);
        end

        % Remove the first cell
        parameter_sim_CellArray(1) = [];

    case 'e_pooling'
        %% 5. Entropy State conditioned
        %%%% e_pooling: (alpha, tau_prior, z_target [can be 'mean','latest' or scalar])%%%%

        % storage array
        parameter_sim_CellArray = cell(1,1);

        Alpha_range = linspace(0.1, 1, num_parameters); % Double (needs to be at least 0.1 to allow convergence in MAXSR)
        Tau_Prior_range = linspace(3, 36, num_parameters); % Double

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Loop through all signals %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for signal = 1:num_signals

            %%% Generate configurations %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % get specific Z_target per signal based on bounds of signal
            [Sig_min,Sig_max] = bounds(signal_series{:,signal});
            Z_target_range = linspace(Sig_min, Sig_max, num_parameters); % Generated based on signal range (range between min and max)

            hyperparameter_set = [Z_target_range; Alpha_range; Tau_Prior_range];

            % Initialize an empty matrix to store configurations
            num_configs = size(hyperparameter_set, 2)^size(hyperparameter_set, 1);
            parameter_configMatrix = zeros(num_configs, size(hyperparameter_set, 1));

            % Nested loops to generate all configurations based on the size of hyperparameter_set
            col = 1;
            for config1 = 1:size(hyperparameter_set, 2)
                for config2 = 1:size(hyperparameter_set, 2)
                    for config3 = 1:size(hyperparameter_set, 2)
                        % Store the current combination of values in the configurations matrix
                        parameter_configMatrix(col, :) = [hyperparameter_set(1, config1); hyperparameter_set(2, config2); hyperparameter_set(3, config3)];
                        col = col + 1;
                    end
                end
            end

            base_backtestObject.Signals = signal_series(:,signal);


            %%%%%%%%%%%%%%% Apply configuration in simulations %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% Create storage for simulation objects
            parameter_num_configs = size(parameter_configMatrix, 1);
            parameter_simulationCellArray = cell(1, parameter_num_configs);

            % Create a cell array to store the test objects
            parameter_test_class1Array = cell(1, parameter_num_configs);

            % Populate the cell array with initial objects
            parfor config = 1:parameter_num_configs
                iterativeClass = copy(base_backtestObject);  % Use the copy method
                iterativeClass.HSFPparameters.Z_target = parameter_configMatrix(config, 1);
                iterativeClass.HSFPparameters.Alpha = parameter_configMatrix(config, 2);
                iterativeClass.HSFPparameters.Tau_prior = parameter_configMatrix(config, 3);

                % assign object
                parameter_test_class1Array{config} = iterativeClass;
            end

            % Use parfor to run backtests
            parfor config = 1:parameter_num_configs
                parameter_simulationCellArray{config} = OOPbacktest_analysis(parameter_test_class1Array{config});
            end

            % Accumulate results each signal loop
            parameter_sim_CellArray = horzcat(parameter_sim_CellArray,parameter_simulationCellArray);
        end

        % Remove the first cell
        parameter_sim_CellArray(1) = [];

    case 'ew_ensemble'
        %% 4. Equally Weighted Time & State conditioned %%%%%%%%%%%%%%%%%%%%%
        %%%% ew_ensemble: (alpha, tau_prior, z_target [can ONLY be 'mean' or 'latest']) %%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Create Signal Combos %%%%%%
        % Assuming you have a cell array 'variables' with your 12 variables
        variables = 1:width(signal_series);
        % Set hard limit of 3 as max number of varaibles in a given combo
        num_variables = 3;

        % Initialize an empty cell array to store combinations
        all_sig_combos = cell(1, 1);

        % Nested loops to generate all combinations
        for num_combinations = 2:num_variables
            combinations = nchoosek(variables, num_combinations);

            % Append the current combinations to the overall list
            all_sig_combos = [all_sig_combos; combinations];
        end

        % Remove the first cell
        all_sig_combos(1) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % storage array
        parameter_sim_CellArray = cell(1,1);


        %%% Generate configurations %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Alpha_range = linspace(0.1, 1, num_parameters); % Double (needs to be at least 0.1 to allow convergence in MAXSR)
        Tau_Prior_range = linspace(3, 36, num_parameters); % Double

        hyperparameter_set = [Alpha_range; Tau_Prior_range];

        % Initialize an empty matrix to store configurations
        num_configs = size(hyperparameter_set, 2)^size(hyperparameter_set, 1);
        parameter_configMatrix = zeros(num_configs, size(hyperparameter_set, 1));

        % Nested loops to generate all configurations based on the size of hyperparameter_set
        col = 1;
        for config1 = 1:size(hyperparameter_set, 2)
            for config2 = 1:size(hyperparameter_set, 2)
                % Store the current combination of values in the configurations matrix
                parameter_configMatrix(col, :) = [hyperparameter_set(1, config1); hyperparameter_set(2, config2)];
                col = col + 1;
            end
        end

        % Account for the binary option for gamma of kernel:
        % z_taget = 1 is 'mean' mu(Z)
        % z_taget = 2 is 'latest' (z_T)
        configurations_mean1 = [parameter_configMatrix, ones(size(parameter_configMatrix, 1), 1)];
        configurations_latest2 = [parameter_configMatrix, 2*ones(size(parameter_configMatrix, 1), 1)];
        parameter_configMatrix = [configurations_mean1; configurations_latest2];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Loop through each possible sig combo %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through the cell array based on num of variables per combo
        for num = 1:num_variables-1
            % Loop over each combo
            for combo = 1:length(all_sig_combos{num,1})
                sig_idxs = all_sig_combos{num,1}(combo,:);
                % Subsetting the timetable
                subsetT = signal_series(:,sig_idxs);

                % update signal per loop
                base_backtestObject.Signals = subsetT;

                %%%%%%%%%%%%%%% Apply configuration in simulations %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%% Create storage for simulation objects
                parameter_num_configs = size(parameter_configMatrix, 1);
                parameter_simulationCellArray = cell(1, parameter_num_configs);

                % Create a cell array to store the test objects
                parameter_test_class1Array = cell(1, parameter_num_configs);

                % Populate the cell array with initial objects
                parfor config = 1:parameter_num_configs
                    iterativeClass = copy(base_backtestObject);  % Use the copy method
                    iterativeClass.HSFPparameters.Alpha = parameter_configMatrix(config, 1);
                    iterativeClass.HSFPparameters.Tau_prior = parameter_configMatrix(config, 2);
                    if parameter_configMatrix(config, 3) == 1
                        iterativeClass.HSFPparameters.Z_target = 'mean';
                    elseif parameter_configMatrix(config, 3) == 2
                        iterativeClass.HSFPparameters.Z_target = 'latest';
                    end

                    % assign object
                    parameter_test_class1Array{config} = iterativeClass;
                end

                % Use parfor to run backtests
                parfor config = 1:parameter_num_configs
                    parameter_simulationCellArray{config} = OOPbacktest_analysis(parameter_test_class1Array{config});
                end

                % Accumulate results each signal loop
                parameter_sim_CellArray = horzcat(parameter_sim_CellArray,parameter_simulationCellArray);

            end


        end

        % Remove the first cell
        parameter_sim_CellArray(1) = [];


    case 'cb_ensemble'
        %% 4. Equally Weighted Time & State conditioned %%%%%%%%%%%%%%%%%%%%%
        %%%% cb_ensemble: (alpha, tau_prior, z_target [can ONLY be 'mean' or 'latest'] ,
        % ensemble_wt_method)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Create Signal Combos %%%%%%
        % Assuming you have a cell array 'variables' with your 12 variables
        variables = 1:width(signal_series);
        % Set hard limit of 3 as max number of varaibles in a given combo
        num_variables = 3;

        % Initialize an empty cell array to store combinations
        all_sig_combos = cell(1, 1);

        % Nested loops to generate all combinations
        for num_combinations = 2:num_variables
            combinations = nchoosek(variables, num_combinations);

            % Append the current combinations to the overall list
            all_sig_combos = [all_sig_combos; combinations];
        end

        % Remove the first cell
        all_sig_combos(1) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % storage array
        parameter_sim_CellArray = cell(1,1);


        %%% Generate configurations %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Alpha_range = linspace(0.1, 1, num_parameters); % Double (needs to be at least 0.1 to allow convergence in MAXSR)
        Tau_Prior_range = linspace(3, 36, num_parameters); % Double

        hyperparameter_set = [Alpha_range; Tau_Prior_range];

        % Initialize an empty matrix to store configurations
        num_configs = size(hyperparameter_set, 2)^size(hyperparameter_set, 1);
        parameter_configMatrix = zeros(num_configs, size(hyperparameter_set, 1));

        % Nested loops to generate all configurations based on the size of hyperparameter_set
        col = 1;
        for config1 = 1:size(hyperparameter_set, 2)
            for config2 = 1:size(hyperparameter_set, 2)
                % Store the current combination of values in the configurations matrix
                parameter_configMatrix(col, :) = [hyperparameter_set(1, config1); hyperparameter_set(2, config2)];
                col = col + 1;
            end
        end

        % Account for the binary option for gamma of kernel:
        % z_taget = 1 is 'mean' mu(Z)
        % z_taget = 2 is 'latest' (z_T)
        configurations_mean1 = [parameter_configMatrix, ones(size(parameter_configMatrix, 1), 1)];
        configurations_latest2 = [parameter_configMatrix, 2*ones(size(parameter_configMatrix, 1), 1)];
        parameter_configMatrix = [configurations_mean1; configurations_latest2];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Loop through each possible sig combo %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through the cell array based on num of variables per combo
        for num = 1:num_variables-1
            % Loop over each combo
            for combo = 1:length(all_sig_combos{num,1})
                sig_idxs = all_sig_combos{num,1}(combo,:);
                % Subsetting the timetable
                subsetT = signal_series(:,sig_idxs);

                % update signal per loop
                base_backtestObject.Signals = subsetT;

                %%%%%%%%%%%%%%% Apply configuration in simulations %%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%% Create storage for simulation objects
                parameter_num_configs = size(parameter_configMatrix, 1);
                parameter_simulationCellArray = cell(1, parameter_num_configs);

                % Create a cell array to store the test objects
                parameter_test_class1Array = cell(1, parameter_num_configs);

                % Populate the cell array with initial objects
                parfor config = 1:parameter_num_configs
                    iterativeClass = copy(base_backtestObject);  % Use the copy method
                    iterativeClass.HSFPparameters.Alpha = parameter_configMatrix(config, 1);
                    iterativeClass.HSFPparameters.Tau_prior = parameter_configMatrix(config, 2);
                    if parameter_configMatrix(config, 3) == 1
                        iterativeClass.HSFPparameters.Z_target = 'mean';
                    elseif parameter_configMatrix(config, 3) == 2
                        iterativeClass.HSFPparameters.Z_target = 'latest';
                    end

                    % assign object
                    parameter_test_class1Array{config} = iterativeClass;
                end

                % Use parfor to run backtests
                parfor config = 1:parameter_num_configs
                    parameter_simulationCellArray{config} = OOPbacktest_analysis(parameter_test_class1Array{config});
                end

                % Accumulate results each signal loop
                parameter_sim_CellArray = horzcat(parameter_sim_CellArray,parameter_simulationCellArray);

            end


        end

        % Remove the first cell
        parameter_sim_CellArray(1) = [];

    otherwise
        % Handle the case when none of the above conditions are met
        error('Invalid method specified.');
end

end




