% quasiDiagSort sub-function of HRP estimate  

function sortedIdx = quasiDiagSort(link)
% quasiDiagSort orders the nodes according to the similarities
% specified by linkage matrix.

% Copyright 2019 The MathWorks, Inc.
numLeafNodes = size(link, 1) + 1;
rootGroupNodeId = 2*numLeafNodes-1;
sortedIdx = getLeafNodesInGroup(rootGroupNodeId, link);
end
