% qdiag sub-function of HRP estimate performs the quasi-diagonalization 

function sortedIdx = qdiag(link)
% qdiag orders the nodes according to the similarities
% specified by linkage matrix.

% Copyright 2019 The MathWorks, Inc.
nLeafNodes = size(link, 1) + 1;
rootNodeId = 2*nLeafNodes-1;
sortedIdx = leafnodes(rootNodeId, link);
end
