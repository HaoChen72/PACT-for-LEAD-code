function results = group_comparison(X, y, groupIdxCell, groupNameCell)
% GROUP_STACKED_TEST
%   Fror each group, all (subject x feature) normalized values are pooled into 
%   one vector per class, then compared via Welch two-sample t-test.
%   Also returns group means and Cohen's d.
%
% Inputs:
%   X             : [N x P] feature matrix
%                   Pre-normalize with: X = normalize(X, 1, 'range');
%   y             : [N x 1] binary labels (0/1)
%   groupIdxCell  : {1 x G} cell, each = vector of column indices
%   groupNameCell : {1 x G} cell of group name strings (optional)

if nargin < 3
    error('Need at least X, y, and groupIdxCell.');
end
if nargin < 4 || isempty(groupNameCell)
    G = numel(groupIdxCell);
    groupNameCell = arrayfun(@(g) sprintf('group_%d',g), 1:G, 'UniformOutput', false);
end

X = normalize(X,1,'range');
y   = y(:);
N   = size(X,1);
G   = numel(groupIdxCell);
ind0 = find(y == 0);
ind1 = find(y == 1);

if numel(y) ~= N
    error('X and y must have the same number of rows.');
end

% Preallocate
nFeat    = nan(G,1);
mean0    = nan(G,1);
mean1    = nan(G,1);
meanDiff = nan(G,1);
cohensD  = nan(G,1);
pRaw     = nan(G,1);
effN0    = nan(G,1);
effN1    = nan(G,1);

for g = 1:G
    idx = groupIdxCell{g}(:)';
    nFeat(g) = numel(idx);

    % Stack: subjects x features -> single column, then remove NaNs
    v0 = X(ind0, idx); v0 = v0(:); v0 = v0(~isnan(v0));
    v1 = X(ind1, idx); v1 = v1(:); v1 = v1(~isnan(v1));

    effN0(g) = numel(v0);
    effN1(g) = numel(v1);

    if numel(v0) < 2 || numel(v1) < 2
        continue;
    end

    mean0(g)    = mean(v0);
    mean1(g)    = mean(v1);
    meanDiff(g) = mean1(g) - mean0(g);

    % Welch t-test (same as default ttest2)
    [~, pRaw(g)] = ttest2(v0, v1, 'Vartype', 'unequal');

    % Cohen's d on stacked vectors (pooled SD)
    n0 = numel(v0); n1 = numel(v1);
    sp = sqrt(((n0-1)*var(v0,0) + (n1-1)*var(v1,0)) / (n0+n1-2));
    if sp > 0
        cohensD(g) = meanDiff(g) / sp;
    else
        cohensD(g) = 0;
    end
end

% BH-FDR correction across groups
pBH = bh_fdr(pRaw);

results = table( ...
    string(groupNameCell(:)), ...
    nFeat(:), ...
    effN0(:), ...
    effN1(:), ...
    mean0(:), ...
    mean1(:), ...
    meanDiff(:), ...
    cohensD(:), ...
    pRaw(:), ...
    pBH(:), ...
    'VariableNames', { ...
        'Group', ...
        'nFeatures', ...
        'EffN_group0', ...
        'EffN_group1', ...
        'Mean_group0', ...
        'Mean_group1', ...
        'MeanDiff_1minus0', ...
        'CohensD', ...
        'p_raw', ...
        'p_BH'} );
end

% -------- Local helper function --------
function q = bh_fdr(p)
p  = p(:);
m  = numel(p);
[ps, idx]  = sort(p);
qsorted    = ps .* m ./ (1:m)';
qsorted    = flipud(cummin(flipud(qsorted)));
qsorted(qsorted > 1) = 1;
q          = nan(m,1);
q(idx)     = qsorted;
end
