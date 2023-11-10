function [corr1, clstrs, silh] = clusterKMeansBase(corr0, max_K, n_init,options)
% This function performs a ** FIRST-pass ** estimate of E[K] for an unsupervised
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
if nargin < 3
    n_init = 10;
end

if nargin < 4
    options = [];
end


%%%%%%
corr0_table = array2table(corr0);

% Calculate distance matrix (dist) and initialize silh
dist = sqrt((1 - fillmissing(corr0, 'constant', 0)) / 2);
silh_vec = []; % Creates an empty vector

for init = 1:n_init
    for i = 3:max_K
        % Set random number generator for consistency:
        rng(1);
        [idx,~,~,~] = kmeans(dist, i,'Distance','sqeuclidean','Options' ...
            ,options,'MaxIter',1000,...
            'Display','off','Replicates',1);

        rng(1);
        silh_ = silhouette(dist, idx);
        % Stores current and previous quality measure for comparision
        stat = [mean(silh_)/std(silh_), mean(silh_vec)/std(silh_vec)];


        if isnan(stat(2)) || stat(1) > stat(2)
            silh_vec = silh_;
          
            n_clusters = length(unique(idx)); % establishes how many clusters kmeans chose
            [~, newIdx] = sort(idx); % shows where the return series lies on a sorted inde listed based on its assigned cluster
            corr1 = corr0(newIdx, newIdx); % The r & c of the input corr matrix reordered based on cluster assignments.

            % Define DICTIONARY
            clstrs = configureDictionary("double","cell");

            % Iterate through each cluster: KEY
            for j = 1:n_clusters

                % Find the indices of variables in the current cluster
                clusterIndices = find(idx == j);

                % Store the variable names in the j-th cell
                clstrs{j} = clusterIndices;
            end
        end
    end
end

silh = silh_vec';
end
