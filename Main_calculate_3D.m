%% Automated 3D Visualization Generation and Vessel Parameter Calculation
% Retained metrics only: BranchNum, BranchPoints, sumDM, sumSOAM, MaxLen, MeanBranchLen
% Author: Dong Shirui (Integrated Modification)

% FF denotes the forefoot region
% MB denotes the combined midfoot and backfoot region

%% Input: Set root path
root_path = 'C:\Users\chnzi\Desktop\PACT-for-LEAD-code-main\exampleData\'; % Modify to your data path
clc; close all force; warning('off','images:initSize:adjustingMag');

%% Path validation and directory tree construction
if ~isfolder(root_path)
    error('Root path does not exist: %s', root_path);
end
primary_list = dir(fullfile(root_path, '*'));
primary_list = primary_list([primary_list.isdir]);
primary_list = primary_list(~ismember({primary_list.name}, {'.','..'}));
total_primary = numel(primary_list);
fprintf('Found %d primary folders\n', total_primary);

%% ===== Summary Table Column Definitions =====
varNames = { ...
    'FolderPath', ...
    'Global_BranchNum','Global_BranchPoints','Global_sum_DM','Global_sum_SOAM','Global_MaxLen','Global_MeanBranchLen', ...
    'FF_Large_BranchNum','FF_Large_BranchPoints','FF_Large_sum_DM','FF_Large_sum_SOAM','FF_Large_MaxLen','FF_Large_MeanBranchLen', ...
    'FF_Small_BranchNum','FF_Small_BranchPoints','FF_Small_sum_DM','FF_Small_sum_SOAM','FF_Small_MaxLen','FF_Small_MeanBranchLen', ...
    'MB_Large_BranchNum','MB_Large_BranchPoints','MB_Large_sum_DM','MB_Large_sum_SOAM','MB_Large_MaxLen','MB_Large_MeanBranchLen', ...
    'MB_Small_BranchNum','MB_Small_BranchPoints','MB_Small_sum_DM','MB_Small_sum_SOAM','MB_Small_MaxLen','MB_Small_MeanBranchLen' ...
    };
varTypes = ["string", repmat("double",1,numel(varNames)-1)];
results_summary = table('Size',[0 numel(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);

%% ===== Main Processing Flow =====
for primary_idx = 1:total_primary
    current_primary = fullfile(primary_list(primary_idx).folder, primary_list(primary_idx).name);
    secondary_path = fullfile(current_primary, 'data');
    if ~isfolder(secondary_path)
        fprintf('[WARNING] Skipping directory missing data: %s\n', current_primary);
        continue;
    end
    tertiary_list = dir(fullfile(secondary_path, '*'));
    tertiary_list = tertiary_list([tertiary_list.isdir]);
    tertiary_list = tertiary_list(~ismember({tertiary_list.name}, {'.','..'}));
    process_tertiary = tertiary_list(1:min(2,end)); % Take first two only (L/R)

    tic;
    for tertiary_idx = 1:numel(process_tertiary)
        %% Clean loop variables but keep necessary external variables
        clearvars -except primary_idx primary_list total_primary tertiary_idx process_tertiary ...
            results_summary root_path current_primary secondary_path tic

        current_tertiary = fullfile(process_tertiary(tertiary_idx).folder, process_tertiary(tertiary_idx).name);
        data_file = fullfile(current_tertiary, 'top_view_seg.mat');
        if ~isfile(data_file)
            fprintf('[ERROR] Data file missing: %s\n', data_file);
            continue;
        end

        try
            temp = load(data_file);
            top_view_seg = temp.top_view_seg;
        catch ME
            fprintf('[LOAD FAILED] %s\nReason: %s\n', data_file, ME.message);
            continue;
        end
        if ~isnumeric(top_view_seg) || ndims(top_view_seg) ~= 2
            warning('top_view_seg is not a 2D matrix, skipping: %s', current_tertiary);
            continue;
        end

        %% ===== Generate Vessel Mask and Skeleton =====
        A = vessel_seg_PA3D_v1(top_view_seg); % Binary vessel mask (1=vessel)
        if ~islogical(A); A = logical(A); end
        ve_skel = bwmorph(A,'skel',Inf);      % Single pixel skeleton
        ve_skel = bwmorph(ve_skel,'spur',5);  % Remove short spurs
        ve_skel_seg = bwareaopen(ve_skel, 9); % Remove skeletons that are too short (<9px)
        [L_skel, vessel_count] = bwlabel(ve_skel_seg, 8); % Connected components on skeleton

        %% ===== Basic Parameters (Whole Foot) =====
        branch_points_map = bwmorph(ve_skel_seg, 'branchpoints');
        total_branch_points = nnz(branch_points_map);
        total_branch_num = total_branch_points + vessel_count;

        ve_props = regionprops(L_skel, 'PixelIdxList','Centroid');

        %% Initialize vessel parameters
        len_k = zeros(1, vessel_count);
        br_pts_k = zeros(1, vessel_count);
        br_num_k = zeros(1, vessel_count);
        dm_k = zeros(1, vessel_count);
        soam_rad_k = zeros(1, vessel_count);        
        branch_len_k = zeros(1, vessel_count);

        %% Calculate vessel parameters
        for k = 1:vessel_count
            % Main Trunk Skeleton Subset
            skel_sub = false(size(A));
            skel_sub(ve_props(k).PixelIdxList) = true;
            total_skel_length = numel(ve_props(k).PixelIdxList);

            % Main Trunk Path and Main Vessel Length
            [mainPath, mainLength] = findLongestPath(skel_sub);  % Requires subfunction findLongestPath
            len_k(k) = mainLength;

            % Branch Calculation
            pix = ve_props(k).PixelIdxList;
            br_pts_k(k) = sum(branch_points_map(pix));
            br_num_k(k) = br_pts_k(k) + 1;

            % Branch Length Calculation
            if br_pts_k(k) > 0
                branch_len_k(k) = max(total_skel_length - mainLength, 0) / br_pts_k(k);
            else
                branch_len_k(k) = 0;
            end

            % Tortuosity DM (Based on main trunk)
            if size(mainPath,1) >= 2
                actual_len = mainLength;
                p_start = mainPath(1,:);
                p_end   = mainPath(end,:);
                eu_dist = sqrt(sum((p_end - p_start).^2));
                if eu_dist > 0
                    dm_k(k) = actual_len / eu_dist - 1;
                else
                    dm_k(k) = 0;
                end
            else
                dm_k(k) = 0;
            end

            % Tortuosity SOAM (Based on main trunk)
            total_angle = 0; valid_points = 0;
            if size(mainPath,1) >= 3
                for j = 2:(size(mainPath,1)-1)
                    v1 = mainPath(j-1,:) - mainPath(j,:);
                    v2 = mainPath(j+1,:) - mainPath(j,:);
                    n1 = norm(v1); n2 = norm(v2);
                    if n1==0 || n2==0, continue; end
                    cth = dot(v1,v2)/(n1*n2);
                    cth = min(max(cth,-1),1);
                    theta = acos(cth);
                    turning = pi - theta;
                    total_angle = total_angle + turning;
                    valid_points = valid_points + 1;
                end
            end
            soam_rad_k(k) = (valid_points>0) * total_angle;
        end

        %% Global Metrics
        sum_curve_dm_global = sum(dm_k, 'omitnan');
        sum_soam_rad_global = sum(soam_rad_k, 'omitnan');
        max_len_global = max(len_k);
        mean_branch_len_global = mean(len_k ./ br_num_k, 'omitnan');

        %% ===== Segmented Calculations =====
        W = size(A,2);
        split_col = max(1, round(W/3));
        seg_col_range = { 1:split_col, (split_col+1):W };

        centroids = cat(1, ve_props.Centroid);
        centroid_cols = round(centroids(:,2));
        centroid_cols = min(max(centroid_cols,1), W);

        is_small = (len_k >= 9) & (len_k <= 99);
        is_large = (len_k >= 100);

        idx_f13 = find(centroid_cols >= seg_col_range{1}(1) & centroid_cols <= seg_col_range{1}(end));
        idx_f13_small = idx_f13(is_small(idx_f13));
        idx_f13_large = idx_f13(is_large(idx_f13));
        idx_l23 = find(centroid_cols >= seg_col_range{2}(1) & centroid_cols <= seg_col_range{2}(end));
        idx_l23_small = idx_l23(is_small(idx_l23));
        idx_l23_large = idx_l23(is_large(idx_l23));

        % Aggregate retained metrics for each segment
        FF_L = aggregate_for_selection(idx_f13_large, len_k, br_num_k, br_pts_k, dm_k, soam_rad_k, branch_len_k);
        FF_S = aggregate_for_selection(idx_f13_small, len_k, br_num_k, br_pts_k, dm_k, soam_rad_k, branch_len_k);
        MB_L = aggregate_for_selection(idx_l23_large, len_k, br_num_k, br_pts_k, dm_k, soam_rad_k, branch_len_k);
        MB_S = aggregate_for_selection(idx_l23_small, len_k, br_num_k, br_pts_k, dm_k, soam_rad_k, branch_len_k);

        %% ===== Write to Summary Table =====
        new_row = { ...
            current_tertiary, ...
            total_branch_num, total_branch_points, sum_curve_dm_global, sum_soam_rad_global, max_len_global, mean_branch_len_global, ...
            FF_L.BranchNum,FF_L.BranchPoints,FF_L.sumDM,FF_L.sumSOAM,FF_L.MaxLen,FF_L.MeanBranchLen, ...
            FF_S.BranchNum,FF_S.BranchPoints,FF_S.sumDM,FF_S.sumSOAM,FF_S.MaxLen,FF_S.MeanBranchLen, ...
            MB_L.BranchNum,MB_L.BranchPoints,MB_L.sumDM,MB_L.sumSOAM,MB_L.MaxLen,MB_L.MeanBranchLen, ...
            MB_S.BranchNum,MB_S.BranchPoints,MB_S.sumDM,MB_S.sumSOAM,MB_S.MaxLen,MB_S.MeanBranchLen ...
            };
        results_summary = [results_summary; new_row];

        clear top_view_seg A ve_skel ve_skel_seg L_skel ve_props ...
            branch_points_map len_k br_pts_k br_num_k dm_k soam_rad_k ...
            centroid_cols branch_len_k
    end
    fprintf('\nProcessed primary folder %d/%d! Time: %.2f seconds\n', primary_idx, total_primary, toc);
end

%% Save Excel Summary
if ~isempty(results_summary)
    excel_file = fullfile(root_path, 'vessel_analysis_3D.xlsx');
    writetable(results_summary, excel_file);
    fprintf('\nResults summary table saved to: %s\n', excel_file);
else
    warning('No result data generated, unable to create Excel summary table');
end

%% ===== Local Subfunction: Aggregate Retained Metrics =====
function out = aggregate_for_selection(idx, len_k, br_num_k, br_pts_k, dm_k, soam_rad_k, branch_len_k)
    if isempty(idx)
        out.BranchNum = 0;
        out.BranchPoints = 0;
        out.sumDM = 0;
        out.sumSOAM = 0;
        out.MaxLen = 0;
        out.MeanBranchLen = 0;
    else
        out.BranchNum = sum(br_num_k(idx));
        out.BranchPoints = sum(br_pts_k(idx));
        out.sumDM = sum(dm_k(idx));
        out.sumSOAM = sum(soam_rad_k(idx));
        out.MaxLen = max(len_k(idx));
        out.MeanBranchLen = mean(branch_len_k(idx), 'omitnan');
    end
end