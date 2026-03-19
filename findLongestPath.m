function [mainPath, mainLength] = findLongestPath(skel)
% Input skel: Subset of a single vessel skeleton (logical matrix)
% Output mainPath: N x 2 pixel coordinates (row, col)
%        mainLength: Main trunk length (pixels)

[H,W] = size(skel);
pix = find(skel);
n_pix = numel(pix);
mainPath = [];
mainLength = 0;

if n_pix < 2
    return;
end

[r_coords, c_coords] = ind2sub([H,W], pix);
r_coords = double(r_coords(:)); c_coords = double(c_coords(:));

% Construct adjacency graph
lin2node = zeros(H*W,1,'int32');
lin2node(pix) = int32(1:n_pix);

% Generate edges
neigh_offsets = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];
ei = zeros(0,1); ej = zeros(0,1); ew = zeros(0,1);
for ii = 1:n_pix
    rr = r_coords(ii); cc = c_coords(ii);
    for oo = 1:8
        r2 = rr + neigh_offsets(oo,1);
        c2 = cc + neigh_offsets(oo,2);
        if r2>=1 && r2<=H && c2>=1 && c2<=W
            lin2 = sub2ind([H,W], r2, c2);
            nid = lin2node(lin2);
            if nid>0 && nid>ii
                ei(end+1,1) = ii; %#ok<AGROW>
                ej(end+1,1) = nid; %#ok<AGROW>
                ew(end+1,1) = sqrt((rr-r2)^2 + (cc-c2)^2); %#ok<AGROW>
            end
        end
    end
end

if isempty(ei)
    return;
end

G = graph(ei, ej, ew, n_pix);

% Find the longest path: Take the maximum value of the distance matrix
D = distances(G);
[mainLength, idx] = max(D(:));
if mainLength==0
    mainPath = [];
    return;
end
[row_idx, col_idx] = ind2sub([n_pix,n_pix], idx);

% Use shortestpath to get the node sequence
nodePath = shortestpath(G, row_idx, col_idx);
mainPath = [r_coords(nodePath), c_coords(nodePath)];
end