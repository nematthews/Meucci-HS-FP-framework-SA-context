function [data_processing_sim_ExcessRet_cellArr] = extractSimExcessRet(sim_cellArr)
% Extract the series of excess returns for both MV and HRP portfolios from 
% each simulated backtest and store them in a single array.
%
%% INPUTS:
% sim_cellArr - array that contains M simulated backtested portfolios where
% each simulation is based on a different configuration of input
% parameters.
% (Type: cell array [1 x M])
%
%% OUTPUT:
% data_processing_sim_ExcessRet_cellArr - cell array with the first row
% consissting of the portfolio type names i.e. "MV" & "HRP" and the second
% row consists of a [T x M] timetable for each portfolio type containing
% the excess returns from all M simulations.
% (Type: cell array [2 x 2])

% Author: Nina Matthews (2023)

% $Revision: 1.1 $ $Date: 2023/12/04 09:28:07 $ $Author: Nina Matthews $

%%%%%%%% MV Simulations %%%
% Initialize an empty T x M timetable
M = length(sim_cellArr);
T = length(sim_cellArr{1}.Realised_tsPRet_TT.Time);
MV_simulation_ExcessRet_TT = timetable(sim_cellArr{1}.Realised_tsPRet_TT(2:end, :).Time);

% Iterate through the cells in the cell array
for m = 1:M
    % Extract the T x 1 double object from the cell
    data = sim_cellArr{m}.ExcessReturns(:,2);
    
    % Convert the T x 1 double to a table and add it as a column to the new timetable
    MV_simulation_ExcessRet_TT = addvars(MV_simulation_ExcessRet_TT, data);
end

MV_simulation_ExcessRet_TT.Properties.Description =  "Excess Returns extracted from MV optimisations run over " + num2str(M) + " trials. Frequency: MONTHLY. Last updated: " + datestr(now, 'yyyy-mm-dd') + ". Intended storage: cell array along with HRP excess returns in intended file name: DATA-SIM-BKTST_MV_HRP-EXCESS_RET-M-CELL_ARR.mat";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%%%%%%% HRP Simulations %%%
HRP_simulation_ExcessRet_TT = timetable(sim_cellArr{1}.Realised_tsPRet_TT(2:end, :).Time);

% Iterate through the cells in the cell array
for m = 1:M
    % Extract the T x 1 double object from the cell
    data = sim_cellArr{m}.ExcessReturns(:,4);
    
    % Convert the T x 1 double to a table and add it as a column to the new timetable
    HRP_simulation_ExcessRet_TT = addvars(HRP_simulation_ExcessRet_TT, data);
end

HRP_simulation_ExcessRet_TT.Properties.Description = "Excess Returns extracted from HRP optimisations run over " + num2str(M) + " trials. Frequency: MONTHLY. Last updated: " + datestr(now, 'yyyy-mm-dd') + ". Intended storage: cell array along with MV excess returns in intended file name: DATA-SIM-BKTST_MV_HRP-EXCESS_RET-M-CELL_ARR.mat";


% Store Realised Returns in single cell array for convenience
data_processing_sim_ExcessRet_cellArr = cell(2,2);
data_processing_sim_ExcessRet_cellArr{1,1} = "MV Sim ExcessRet TT";
data_processing_sim_ExcessRet_cellArr{2,1} = MV_simulation_ExcessRet_TT;
data_processing_sim_ExcessRet_cellArr{1,2} = "HRP Sim ExcessRet TT";
data_processing_sim_ExcessRet_cellArr{2,2} = HRP_simulation_ExcessRet_TT;

end