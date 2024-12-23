function ticker_timetable = importfields2keep(fields,fileName,data_start_cell)
% Loops through sheets of excel doc, where each sheet stores
% individual ticker data to import selected fields from all sheets
%
%% INPUTS:
% fields = required fields to store eg. 'PX_LAST' or 'CHG_PCT_1M'
% (type: cell array)
%
% fileName = file and its directory path
% (type: string)
%
% data_start_cell =  Excel file cell reference that maps to the first
% instance of the data that you want to import. Eg. if t he first date sits
% in cell A7, you would specificy "A7".
% (type: string)



% Author: Nina Matthews (2023)

% $Revision: 1.2 $ $Date: 2023/05/04 13:11:05 $ $Author: Nina Matthews $

%% 

% define columns to keep from each sheet
fieldsToKeep = fields;

% create sheetName var to loop through in the file
excelSheetNames = sheetnames(fileName);
% Preallocate cell array for storage before looping
ticker_timetable{numel(excelSheetNames)} = [];
% Load the dataset by sheet
for sheet = 1:numel(excelSheetNames)
    ticker_timetable{sheet} = readtimetable(fileName,'Sheet',excelSheetNames{sheet}, ...
        VariableNamingRule='preserve', VariableNamesRange = data_start_cell, ...
        VariableDescriptionsRange = 'G4',Range=data_start_cell);
    % Assign the ticker as the timetable Description
    ticker_timetable{sheet}.Properties.Description = excelSheetNames{sheet};
    % from all columns remove all but 'varsToKeep'
    fieldsToRemove = setdiff(ticker_timetable{sheet}.Properties.VariableNames, fieldsToKeep);
    if ~isempty(fieldsToRemove)
        ticker_timetable{sheet} = removevars(ticker_timetable{sheet}, fieldsToRemove);
    end
end