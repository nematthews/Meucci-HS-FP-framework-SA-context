function w = hrpbisect(cov, sortedIdx)
% hrpbisect starts with the sorted list from the hierarchical
% clustering tree, and recursively runs bisect allocation on each
% sublist. 


len = numel(sortedIdx);  % can be less than the number of assets
nAssets = size(cov, 1);
w = ones(nAssets, 1);

% base case
if len <= 2
    w(sortedIdx) = ivp(cov(sortedIdx, sortedIdx));
    return;
end

% below demonstrates the bottom-up and top-down approach
% split the list into two:
mid = floor(len/2);
left = sortedIdx(1:mid);
right = sortedIdx(mid+1:len);

% 1) bottom up to find the variance of the two groups
wgtL = ivp(cov(left, left));
wgtR = ivp(cov(right, right));

varLeft = wgtL'*cov(left, left)*wgtL;
varRight = wgtR'*cov(right, right)*wgtR;

% 2) top down rescaling of asset weights in inverse proportion to the two groups' variances
alpha= 1-varLeft/(varLeft+varRight);
w(left) = alpha*w(left);
w(right) = (1-alpha)*w(right);

% recursively run this for the two sublists
w = w.*hrpbisect(cov, left);
w = w.*hrpbisect(cov, right);
end