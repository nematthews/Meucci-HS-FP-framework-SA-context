function [corrNew, clstrsNew, silhNew] = makeNewOutputs(corr0, clstrs, clstrs2)
% 
%
%% INPUTS:
%
% corr0 - correlation matrix of M simulated trails' realised returns.
% (type: double, [M x M])
%
% clstrs - 
%
% clstrs2 - 
%%
    clstrsNew = containers.Map;
    newIdx = [];

    % Populate clstrsNew with clusters from clstrs
    for i = keys(clstrs)
        clstrsNew(clstrsNew.Count + 1) = clstrs(i);
    end

    % Populate clstrsNew with clusters from clstrs2
    for i = keys(clstrs2)
        clstrsNew(clstrsNew.Count + 1) = clstrs2(i);
    end

    % Concatenate values in clstrsNew and update newIdx
    concatenatedValues = vertcat(values(clstrsNew).');

    % Map newIdx to the concatenated values
    for i = 1:length(concatenatedValues)
        newIdx = [newIdx, concatenatedValues{i}];
    end

    % Extract new correlation matrix
    corrNew = corr0(newIdx, newIdx);

    % Calculate the distance matrix
    dist = sqrt((1 - fillmissing(corr0, 'constant', 0)) / 2);

    % Initialize kmeans_labels
    kmeans_labels = zeros(1, size(dist, 2));

    % Assign labels to data points based on clusters in clstrsNew
    for i = keys(clstrsNew)
        idxs = cellfun(@(k) find(strcmp(dist.Properties.RowNames, k)), clstrsNew(i));
        kmeans_labels(idxs) = i;
    end

    % Calculate silhouette values
    silhNew = silhouette(dist, kmeans_labels);

    % Create a MATLAB table for silhNew with the correct row names
    silhNew = array2table(silhNew.', 'VariableNames', dist.Properties.RowNames);

    % Convert the table to a containers.Map
    silhNew = containers.Map(silhNew.Properties.VariableNames, table2cell(silhNew));

    % Return the results
    corrNew = corrNew{:, :};
    clstrsNew = clstrsNew.Values;
    silhNew = silhNew;
end
