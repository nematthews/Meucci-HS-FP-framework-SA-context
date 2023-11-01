function [corr1, clstrs, silh] = clusterKMeansBase(corr0, max_K, n_init,options)
% This function performs a ** first-pass ** estimate of E[K] for an unsupervised
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
% corr0 - correlation matrix of M simulated trails' realised returns
%
% max_K - user defined maximum number of clusters to generate by K-mean.
%
% n_init - number of times the clustering process needs to be repeated.
% Type

%% 
    % Check if maxNumClusters and n_init are provided, otherwise use default values
    if nargin < 2
        max_K = 10;
    end
    if nargin < 3
        n_init = 10;
    end

    corr0_table = array2table(corr0);
    
    % Calculate distance matrix (dist) and initialize silh
    dist = sqrt((1 - fillmissing(corr0, 'constant', 0)) / 2);
    silh_vec = []; % Creates an empty vector

for init = 1:n_init
        for i = 3:max_K
           
           [idx,~,~,~] = kmeans(dist, i,'Distance','sqeuclidean','Options' ...
                            ,options,'MaxIter',100,...
                            'Display','final','Replicates',1);

            silh_ = silhouette(dist, idx);
            % Stores current and previous quality measure for comparision
            stat = [mean(silh_)/std(silh_), mean(silh_vec)/std(silh_vec)];

            
            if isnan(stat(2)) || stat(1) > stat(2)
                silh_vec = silh_;
                % [idx_best,C_best,sumd_best,D_best] = [idx,C,sumd,D];
                n_clusters = length(unique(idx)); % establishes how many clusters kmeans chose
                [~, newIdx] = sort(idx); % shows where the return series lies on a sorted inde listed based on its assigned cluster
                corr1 = corr0(newIdx, newIdx); % The r & c of the input corr matrix reordered based on cluster assignments.
                
               % Define DICTIONARY
                clstrs = dictionary;

                % Iterate through each cluster: KEY
                for j = 1:n_clusters

                    % Find the indices of variables in the current cluster
                    clusterIndices = find(idx == j);
                    
                    % Get the corresponding variable names: VALUES
                    corr0_table.Properties.VariableNames = compose(string(1:width(corr0_table)));
                    clusterVariableNames = corr0_table.Properties.VariableNames(clusterIndices);
                    
                    % Store the variable names in the j-th cell
                    clstrs{j} = clusterVariableNames;
                end
            end
        end
end
    
    silh = silh_vec';
end
