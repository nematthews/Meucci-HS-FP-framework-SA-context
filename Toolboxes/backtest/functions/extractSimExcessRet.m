function [data_processing_sim_ExcessRet_cellArr] = extractSimExcessRet(sim_cellArr)
% Calculate the geometric average of input data, allows for annualisation
% or adjusting for other time periods using f. Typically used to calculate
% portfolio returns for annualising, if so "f" = 12 if returns are monthly. 
%
%
%% INPUTS:
% sim_cellArr - 
% (Type: timetable [T x J])
%
% f - (Optional) number of periods within a yr (e.g 12 for monthly, 
% 4 for quarterly) - Default = 1 i.e not annualised returns
% (type: scalar)  

% Author: Nina Matthews (2023)

% $Revision: 1.1 $ $Date: 2023/013/09 14:18:23 $ $Author: Nina Matthews $

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