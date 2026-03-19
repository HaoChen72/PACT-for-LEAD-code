function rowParams = assembleRowParams(d1,d2,d3,s1,s2,s3,a1,a2,a3,vesselStats, folderName)
% assembleRowParams (Revised Version)
% Output: Original segment-level 60 columns + Whole-foot 16 columns 
% (Correct logic for summation/density/weighting)

rowParams = table();
rowParams.FolderName = folderName;

% Segment-level arrays
densities = [d1, d2, d3];
areas_seg  = [a1, a2, a3];
strengths_seg = [s1, s2, s3];

% --- Segment-level Output (Maintaining original logic) ---
for s = 1:3
    seg = vesselStats(s);

    % Metric extraction
    small_n = seg.small_num;
    med_n   = seg.medium_num;
    large_n = seg.large_num;
    total_n = seg.total_num;

    small_area_per = seg.small_area;   
    med_area_per   = seg.medium_area;
    large_area_per = seg.large_area;

    small_strength = seg.small_strength; 
    med_strength   = seg.medium_strength;
    large_strength = seg.large_strength;

    small_width = seg.small_width;
    med_width   = seg.medium_width;
    large_width = seg.large_width;

    % Aggregate classes
    area_classes = [small_area_per, med_area_per, large_area_per];
    strength_classes = [small_strength, med_strength, large_strength];
    class_nums = [small_n, med_n, large_n];
    class_mean_width = [small_width, med_width, large_width];

    total_area = areas_seg(s);
    total_density = densities(s);
    if total_area > 0
        density_classes = area_classes ./ total_area * total_density;
    else
        density_classes = nan(1,3);
    end

    % Output small/medium/large
    types = {'small','medium','large'};
    for t = 1:3
        prefix = ['seg', num2str(s), '_', types{t}];
        rowParams.([prefix,'_density'])  = density_classes(t);
        rowParams.([prefix,'_strength']) = strength_classes(t);
        rowParams.([prefix,'_area'])      = area_classes(t);
        rowParams.([prefix,'_num'])       = class_nums(t);
        rowParams.([prefix,'_width'])     = class_mean_width(t);
    end

    % total class (within segment)
    prefix = ['seg',num2str(s),'_total'];
    total_area_out = total_area;
    total_density_out = total_density;
    total_strength_out = strengths_seg(s);
    total_num_out = total_n;
    if sum(area_classes) > 0
        total_width_out = sum(class_mean_width .* area_classes) / sum(area_classes);
    else
        % If area info is missing, fall back to num-weighted or simple average
        if sum(class_nums) > 0
            total_width_out = sum(class_mean_width .* class_nums) / sum(class_nums);
        else
            total_width_out = mean(class_mean_width, 'omitnan');
        end
    end

    rowParams.([prefix,'_density'])  = total_density_out;
    rowParams.([prefix,'_strength']) = total_strength_out;
    rowParams.([prefix,'_area'])      = total_area_out;
    rowParams.([prefix,'_num'])       = total_num_out;
    rowParams.([prefix,'_width'])     = total_width_out;
end

% Target: Output density, area, num, width for small/medium/large/total (16 columns total)
% Rules:
%   area, num: Sum of the three segments
%   total density = total_vessel_pixels / total_foot_pixels
%   class density = class_vessel_pixels / total_foot_pixels
%   total_foot_pixels: Recovered from each segment's (area_seg / density_seg)
%   width: Area-weighted average for the same class (fallback to num-weighted or simple average)

% 1) foot_pixel per segment (may contain NaN/0)
foot_pixels_per_seg = nan(1,3);
for s = 1:3
    if ~isnan(densities(s)) && densities(s) ~= 0
        foot_pixels_per_seg(s) = areas_seg(s) / densities(s); % derived from density = area / foot_pixels
    else
        foot_pixels_per_seg(s) = NaN;
    end
end

% 2) Whole-foot total vessel area (pixels) and total foot pixels (estimated)
total_area_all = sum(areas_seg, 'omitnan');
total_foot_pixels_all = sum(foot_pixels_per_seg, 'omitnan');

% If foot_pixels cannot be recovered (all segments are NaN), fall back to NaN
if isempty(total_foot_pixels_all) || total_foot_pixels_all == 0
    total_density_all = NaN;
else
    total_density_all = total_area_all / total_foot_pixels_all;
end

% 3) Total area and total number for each class (sum across three segments)
small_area_all = sum([vesselStats.small_area], 'omitnan');
med_area_all   = sum([vesselStats.medium_area], 'omitnan');
large_area_all = sum([vesselStats.large_area], 'omitnan');

small_num_all = sum([vesselStats.small_num], 'omitnan');
med_num_all   = sum([vesselStats.medium_num], 'omitnan');
large_num_all = sum([vesselStats.large_num], 'omitnan');
total_num_all = sum([vesselStats.total_num], 'omitnan');

% 4) Class density: class_area / total_foot_pixels_all
if ~isnan(total_foot_pixels_all) && total_foot_pixels_all > 0
    small_density_all = small_area_all / total_foot_pixels_all;
    med_density_all   = med_area_all   / total_foot_pixels_all;
    large_density_all = large_area_all / total_foot_pixels_all;
else
    small_density_all = NaN;
    med_density_all   = NaN;
    large_density_all = NaN;
end

% 5) Width: Area-weighted average per class across segments
% Sum (segment class area * segment class width) / Total class area
% If total class area is 0, fall back to num-weighted, then simple average
small_width_num = 0; small_width_den = 0;
med_width_num = 0;   med_width_den = 0;
large_width_num = 0; large_width_den = 0;
for s = 1:3
    seg = vesselStats(s);
    % small
    if ~isnan(seg.small_area) && seg.small_area > 0 && ~isnan(seg.small_width)
        small_width_num = small_width_num + seg.small_width * seg.small_area;
        small_width_den = small_width_den + seg.small_area;
    end
    % med
    if ~isnan(seg.medium_area) && seg.medium_area > 0 && ~isnan(seg.medium_width)
        med_width_num = med_width_num + seg.medium_width * seg.medium_area;
        med_width_den = med_width_den + seg.medium_area;
    end
    % large
    if ~isnan(seg.large_area) && seg.large_area > 0 && ~isnan(seg.large_width)
        large_width_num = large_width_num + seg.large_width * seg.large_area;
        large_width_den = large_width_den + seg.large_area;
    end
end

if small_width_den > 0
    small_width_all = small_width_num / small_width_den;
elseif sum([vesselStats.small_num], 'omitnan') > 0
    % Weighted by number
    small_width_all = sum([vesselStats.small_width] .* [vesselStats.small_num], 'omitnan') / sum([vesselStats.small_num], 'omitnan');
else
    small_width_all = mean([vesselStats.small_width], 'omitnan');
end

if med_width_den > 0
    med_width_all = med_width_num / med_width_den;
elseif sum([vesselStats.medium_num], 'omitnan') > 0
    med_width_all = sum([vesselStats.medium_width] .* [vesselStats.medium_num], 'omitnan') / sum([vesselStats.medium_num], 'omitnan');
else
    med_width_all = mean([vesselStats.medium_width], 'omitnan');
end

if large_width_den > 0
    large_width_all = large_width_num / large_width_den;
elseif sum([vesselStats.large_num], 'omitnan') > 0
    large_width_all = sum([vesselStats.large_width] .* [vesselStats.large_num], 'omitnan') / sum([vesselStats.large_num], 'omitnan');
else
    large_width_all = mean([vesselStats.large_width], 'omitnan');
end

% total width: Weighted by class area (fallback to num-weighted or simple average if total area is 0)
if total_area_all > 0
    total_width_all = (small_width_all * small_area_all + med_width_all * med_area_all + large_width_all * large_area_all) / (small_area_all + med_area_all + large_area_all);
elseif (small_num_all + med_num_all + large_num_all) > 0
    total_width_all = (small_width_all * small_num_all + med_width_all * med_num_all + large_width_all * large_num_all) / (small_num_all + med_num_all + large_num_all);
else
    total_width_all = mean([small_width_all, med_width_all, large_width_all], 'omitnan');
end

% ---- Output 16 Whole-Foot columns (Order: total, large, medium, small; density/area/num/width per class) ----
rowParams.whole_total_density = total_density_all;
rowParams.whole_total_area    = total_area_all;
rowParams.whole_total_num      = total_num_all;
rowParams.whole_total_width   = total_width_all;

rowParams.whole_large_density = large_density_all;
rowParams.whole_large_area    = large_area_all;
rowParams.whole_large_num      = large_num_all;
rowParams.whole_large_width   = large_width_all;

rowParams.whole_medium_density = med_density_all;
rowParams.whole_medium_area    = med_area_all;
rowParams.whole_medium_num      = med_num_all;
rowParams.whole_medium_width   = med_width_all;

rowParams.whole_small_density = small_density_all;
rowParams.whole_small_area    = small_area_all;
rowParams.whole_small_num      = small_num_all;
rowParams.whole_small_width   = small_width_all;

end