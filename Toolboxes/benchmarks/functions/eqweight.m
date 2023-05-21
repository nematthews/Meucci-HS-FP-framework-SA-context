function new_weights = eqweight(asset_ob)
% Equal-weighted portfolio allocation based on number of assets in returns

% If asset_ob is a timetable: Number of assets
if istimetable(asset_ob)
    nAssets = width(asset_ob);
% If asset_ob is an array or matrix: Number of assets
else
    nAssets = size(asset_ob, 2);
end

% Allocate weight of 1 to each asset
new_weights = ones(1, nAssets);

% Divide by the total number of assets to get equal weights
new_weights = new_weights / sum(new_weights);

end