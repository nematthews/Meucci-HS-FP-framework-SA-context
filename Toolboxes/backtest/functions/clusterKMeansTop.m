function [corr1, clstrs, silh] = clusterKMeansTop(corr0, max_K, n_init, options)
% This function performs a ** SECOND-pass ** estimate of E[K] for an unsupervised
% learning approach proposed by Marcos López de Prado & Michael J. Lewis
% (2019).
% The purpose of estimating E[K] is to establish the number of independent
% trails that exist within M simulated trails of backtesting. This function
% should be used in conjunction with simulate_backtest.m that can be used
% to generate M trails given M configurations of J varying
% hyper-parameters.
%
% "For each target number of clusters, we perform a stochastic optimization,
% repeating the clustering operation n_init times. Among all the clustering
% alternatives, we choose the solution that achieves the highest quality
% score, defined as the t-value of the silhouette scores."
%
%% INPUTS:
%
% corr0 - correlation matrix of M simulated trails' realised returns.
% (type: double, [M x M])
%
% max_K - user defined maximum number of clusters to generate by K-mean.
% (type: double)
%
% n_init - number of times the clustering process needs to be repeated.
% (type: double)
%
% options - user defined options that are provided to the kmeans algorithm.
% Used to implement parallel computing abilities.
% E.g.:
% stream = RandStream('mlfg6331_64');  % Random number stream
% options = statset('UseParallel',1,'UseSubstreams',1,...
%     'Streams',stream);
% (type: structure)
%%

% Check if maxNumClusters and n_init are provided, otherwise use default values
if nargin < 2
    max_K = 10;
end
%%%% NOTE: max_K can't exceed number of paths
if max_K > size(corr0,1)
    max_K = size(corr0,1)-1;
end
if nargin < 3
    n_init = 10;
end

if nargin < 4
    options = [];
end

disp("CALL MAIN")
%% Perform base cluster to initialise E[K]
[corr1, clstrs, silh] = clusterKMeansBase(corr0, max_K, n_init,options);
% Calc Tstats per cluster
clusterTstats = dictionary;

for i = 1:length(keys(clstrs))
    % Turn each string element in cells to a double
    series_indx = cellfun(@(x)str2double(x), clstrs{i});
    % Calculate each cluster's Tstat
    format long
    clusterTstats{i} = mean(silh(series_indx)) / std(silh(series_indx));
end

% calc Average tstat across clusters
tStatMean = mean(cell2mat(values(clusterTstats)));

%% extract cluster # that have tstat less than average
redoClusters = find(cell2mat(values(clusterTstats)) < tStatMean);


%% Recluster clusters with low tstats
if length(redoClusters) <= 2
    corrNew = corr1;
    clstrsNew = clstrs;
    silhNew = silh;
    return;
else
    disp("Else entered")
    % Extract all idx values from redo clusters
    keysRedo = [];
    for i = redoClusters
        keysRedo = [keysRedo, clstrs{i}];
    end
    % convert char to double
    keysRedo = cellfun(@(x)str2double(x), keysRedo);

    % Index out the correlations in redo clusters to recluster
    corrTmp = corr0(keysRedo, keysRedo);
    % meanRedoTstat = mean(cell2mat(clusterTstats(redoClusters)));
    disp("Before recursive call")
    %% Recursively call itself

    %%%% NOTE: max_K can't exceed number of paths
    if max_K > size(corrTmp,1)
        max_K = size(corrTmp,1)-1;
    end

    [~, clstrs2, ~] = clusterKMeansTop(corrTmp, max_K, n_init);
    disp("End recursive call")
    % When the call eventually stops as if statement is met:
    % Make new outputs, if necessary

    % Takes all culsters in cluster that are not in RedoClusters and
    % combines it with the new clusters2
    % Iterate over keys of clstrs
    clstrsToKeep = configureDictionary("double","cell");
    for i = keys(clstrs)
        % Check if the key is not in redoClusters
        if ~ismember(i, redoClusters)
            % Add key-value pair to clstrsNew
            clstrsToKeep(i) = clstrs{i};
        end
    end
    
    [corrNew, clstrsNew, silhNew] = makeNewOutputs(corr0, clstrsToKeep, clstrs2);


    newTstatMean = mean(cellfun(@(i) mean(silhNew(clstrsNew(i)))/std(silhNew(clstrsNew(i)), keys(clstrsNew))));

    if newTstatMean <= tStatMean
        dis("If statement entered")
        return;
    else
        corr1 = corrNew;
        clstrs = clstrsNew;
        silh = silhNew;
        dis("Else ENDED")
    end
end
end
