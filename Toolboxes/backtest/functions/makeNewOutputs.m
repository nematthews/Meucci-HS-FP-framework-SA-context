function [corrNew, clstrsNew, silhNew] = makeNewOutputs(corr0, clstrs, clstrs2)
% 
%
%% INPUTS:
%
% corr0 - original correlation matrix of M simulated trails' realised returns.
% (type: double, [M x M])
%
% clstrs - user defined maximum number of clusters to generate by K-mean.
%
% clstrs2 - number of times the clustering process needs to be repeated.
% Type

%% 

    % Initialize variables
    clstrsNew = dictionary;
    newIdx = [];

    % Copy clusters from clstrs
    for i = keys(clstrs)
        clstrsNew(int32(length(clstrsNew.keys) + 1)) = clstrs(i);
    end

    % Copy clusters from clstrs2
    for i = keys(clstrs2)
        clstrsNew(int32(length(clstrsNew.keys) + 1)) = clstrs2(i);
    end

    % Extend newIdx
    cellfun(@(x) newIdx.extend(x), values(clstrsNew));

    % Reorder corr0
    corrNew = corr0(newIdx, newIdx);

    % Calculate distance matrix
    dist = sqrt((1 - fillmissing(corr0, 'constant', 0)) / 2);

    % Initialize kmeans_labels
    kmeans_labels = zeros(1, length(dist.Properties.VariableNames));

    % Assign cluster labels to indices
    for i = keys(clstrsNew)
        idxs = cellfun(@(x) find(ismember(dist.Properties.RowNames, x)), values(clstrsNew(i)));
        kmeans_labels(idxs) = i;
    end

    % Calculate silhouette values
    silhNew = silhouette(dist{:,:}, kmeans_labels);

    % Convert silhouette values to a Series
    silhNew = array2table(silhNew, 'RowNames', dist.Properties.RowNames);
end
