function idxInGroup = leafnodes(NodeId, link)
% leafnodesid finds all leaf nodes for a given group node id
% in a linkage matrix.

% Copyright 2019 The MathWorks, Inc.

N= size(link, 1)+1;
if NodeId>N
    gNodeIds = link(NodeId-N, 1:2);
    idxInGroup = [leafnodes(gNodeIds(1), link), ...
        leafnodes(gNodeIds(2), link)];
else
    idxInGroup = NodeId;
end
end