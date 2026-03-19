%% Automated 3D Visualization Generation and Vessel Parameter Calculation (Segmentation + FWHM + New Basic Parameters)
% Author: Dong Shirui (Integrated Modification)

% F13 denotes the forefoot region
% F23 denotes the combined midfoot and backfoot region

%% Input: Set root path
root_path = 'F:\footData\add2\'; % Modify to your data path
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
    'Global_VesselDensity','Global_MeanIntensity','Global_VesselArea','Global_VesselCount','Global_WidthFWHM','Global_BranchNum','Global_BranchPoints','Global_CurveDM','Global_SOAM_Rad','Global_sum_DM','Global_sum_SOAM','Global_MaxLen','Global_MeanBranchLen', ...
    'F13_Large_Density','F13_Large_Intensity','F13_Large_Area','F13_Large_Count','F13_Large_WidthFWHM','F13_Large_BranchNum','F13_Large_BranchPoints','F13_Large_CurveDM','F13_Large_SOAM_Rad','F13_Large_sum_DM','F13_Large_sum_SOAM','F13_Large_MaxLen','F13_Large_MeanBranchLen', ...
    'F13_Small_Density','F13_Small_Intensity','F13_Small_Area','F13_Small_Count','F13_Small_WidthFWHM','F13_Small_BranchNum','F13_Small_BranchPoints','F13_Small_CurveDM','F13_Small_SOAM_Rad','F13_Small_sum_DM','F13_Small_sum_SOAM','F13_Small_MaxLen','F13_Small_MeanBranchLen', ...
    'L23_Large_Density','L23_Large_Intensity','L23_Large_Area','L23_Large_Count','L23_Large_WidthFWHM','L23_Large_BranchNum','L23_Large_BranchPoints','L23_Large_CurveDM','L23_Large_SOAM_Rad','L23_Large_sum_DM','L23_Large_sum_SOAM','L23_Large_MaxLen','L23_Large_MeanBranchLen', ...
    'L23_Small_Density','L23_Small_Intensity','L23_Small_Area','L23_Small_Count','L23_Small_WidthFWHM','L23_Small_BranchNum','L23_Small_BranchPoints','L23_Small_CurveDM','L23_Small_SOAM_Rad','L23_Small_sum_DM','L23_Small_sum_SOAM','L23_Small_MaxLen','L23_Small_MeanBranchLen' ...
    };
varTypes = ["string", repmat("double",1,numel(varNames)-1)];
results_summary = table('Size',[0 numel(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);

%% ===== Main Processing Flow =====
for primary_idx = 1:total_primary
    current_primary = fullfile(primary_list(primary_idx).folder, primary_list(primary_idx).name);
    secondary_path = fullfile(current_primary, '3D_str');
    if ~isfolder(secondary_path)
        fprintf('[WARNING] Skipping directory missing 3D_str: %s\n', current_primary);
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
        ve_skel = bwmorph(A,'skel',Inf); % Single pixel skeleton
        ve_skel = bwmorph(ve_skel,'spur',5); % Remove short spurs
        ve_skel_seg = bwareaopen(ve_skel, 9); % Remove skeletons that are too short (<9px)
        [L_skel, vessel_count] = bwlabel(ve_skel_seg, 8); % Connected components on skeleton
        [L_area, ~] = bwlabel(A, 8); % Total tissue pixels (non-zero)
        tissue_area_total = nnz(top_view_seg);

        %% ===== Basic Parameters (Whole Foot) =====
        vessel_area_total = nnz(A);
        branch_points_map = bwmorph(ve_skel_seg, 'branchpoints');
        total_branch_points = nnz(branch_points_map);
        total_branch_num = total_branch_points + vessel_count;

        ve_props = regionprops(L_skel, 'Area','PixelList','PixelIdxList','Centroid');
        area_props = regionprops(L_area, 'Area','PixelIdxList','Centroid');

        %% Map skeleton components to area components
        skel2area = zeros(1, vessel_count);
        for k = 1:vessel_count
            pix = ve_props(k).PixelIdxList;
            area_labels = L_area(pix);
            area_labels = area_labels(area_labels>0);
            if isempty(area_labels)
                skel2area(k) = 0;
            else
                skel2area(k) = mode(area_labels);
            end
        end
        %% Initialize vessel parameters
        len_k = zeros(1, vessel_count);
        area_k = zeros(1, vessel_count);
        br_pts_k = zeros(1, vessel_count);
        br_num_k = zeros(1, vessel_count);
        fwhm_k = nan(1, vessel_count);
        dm_k = zeros(1, vessel_count);
        soam_rad_k = zeros(1, vessel_count);        
        branch_len_k = zeros(1, vessel_count);
        [H,W] = size(A);
        top_vec = double(top_view_seg(:));
        %% Calculate vessel parameters
        for k = 1:vessel_count
            % ===== Main Trunk Skeleton Subset =====
            skel_sub = false(size(A));
            skel_sub(ve_props(k).PixelIdxList) = true;
            total_skel_length = numel(ve_props(k).PixelIdxList);  % Total skeleton length

            % ===== Main Trunk Path and Main Vessel Length =====
            [mainPath, mainLength] = findLongestPath(skel_sub);  % Requires subfunction findLongestPath
            len_k(k) = mainLength;  % Replaces original ve_props(k).Area

            % ===== Area Calculation (Keep original logic) =====
            if skel2area(k) > 0
                area_k(k) = numel(area_props(skel2area(k)).PixelIdxList);
            else
                area_k(k) = 0;
            end

            % ===== Branch Calculation (Keep original logic) =====
            pix = ve_props(k).PixelIdxList;
            br_pts_k(k) = sum(branch_points_map(pix));
            br_num_k(k) = br_pts_k(k) + 1;

            % ===== Branch Length Calculation (Optimization) =====
            if br_pts_k(k) > 0
                branch_len_k(k) = max(total_skel_length - mainLength, 0) / br_pts_k(k);
            else
                branch_len_k(k) = 0;  % or NaN
            end

            % ===== Tortuosity DM (Based on main trunk) =====
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

            % ===== Tortuosity SOAM (Based on main trunk) =====
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

            % Radius FWHM
            widths = [];
            sample_step = 10;
            for s = 1:sample_step:size(mainPath,1)
                pt = mainPath(s,:);
                % Local skeleton direction
                if s==1
                    tangent = mainPath(s+1,:) - pt;
                elseif s==size(mainPath,1)
                    tangent = pt - mainPath(s-1,:);
                else
                    tangent = mainPath(s+1,:) - mainPath(s-1,:);
                end
                tangent = tangent / norm(tangent);
                % Normal direction
                normal = [-tangent(2), tangent(1)];
                % Scan pixels along both sides of the normal
                max_dist = 40; % Max scan distance (pixels)
                vals_along = [];
                for d = -max_dist:max_dist
                    r = round(pt(1) + d*normal(1));
                    c = round(pt(2) + d*normal(2));
                    if r>=1 && r<=size(A,1) && c>=1 && c<=size(A,2)
                        vals_along(end+1) = top_view_seg(r,c); 
                    end
                end
                % Call calculateFWHM_v1 to calculate current point width
                if numel(vals_along) >= 4
                    widths(end+1) = calculateFWHM_v1(vals_along); 
                end
            end
            % Main trunk average FWHM
            fwhm_k(k) = mean(widths,'omitnan');
        end

        %% Global Mean Values
        mean_curve_dm_global = mean(dm_k, 'omitnan');
        mean_soam_rad_global = mean(soam_rad_k, 'omitnan');
        sum_curve_dm_global = sum(dm_k, 'omitnan');
        sum_soam_rad_global = sum(soam_rad_k, 'omitnan');
        max_len_global = max(len_k);
        mean_branch_len_global = mean(len_k ./ br_num_k, 'omitnan');
        mean_width_global = mean(fwhm_k, 'omitnan');
        vessel_density_global = vessel_area_total / tissue_area_total;
        mean_intensity_global = mean(nonzeros(double(abs(top_view_seg)).*double(A)));

        %% ===== Segmented Calculations =====
        W = size(A,2);                % Image width
        split_col = max(1, round(W/3));
        seg_col_range = { 1:split_col, (split_col+1):W };

        centroids = cat(1, ve_props.Centroid);  % [row, col]
        centroid_cols = round(centroids(:,2));  % Get column coordinates
        centroid_cols = min(max(centroid_cols,1), W);

        is_small = (len_k >= 9) & (len_k <= 99);
        is_large = (len_k >= 100);

        idx_f13 = find(centroid_cols >= seg_col_range{1}(1) & centroid_cols <= seg_col_range{1}(end));
        idx_f13_small = idx_f13(is_small(idx_f13));
        idx_f13_large = idx_f13(is_large(idx_f13));
        idx_l23 = find(centroid_cols >= seg_col_range{2}(1) & centroid_cols <= seg_col_range{2}(end));
        idx_l23_small = idx_l23(is_small(idx_l23));
        idx_l23_large = idx_l23(is_large(idx_l23));

        % Call subfunction
        F13_L = aggregate_for_selection(idx_f13_large, seg_col_range{1}, A, top_view_seg, skel2area, area_props, len_k, br_num_k, br_pts_k, fwhm_k, dm_k, soam_rad_k, branch_len_k);
        F13_S = aggregate_for_selection(idx_f13_small, seg_col_range{1}, A, top_view_seg, skel2area, area_props, len_k, br_num_k, br_pts_k, fwhm_k, dm_k, soam_rad_k, branch_len_k);
        L23_L = aggregate_for_selection(idx_l23_large, seg_col_range{2}, A, top_view_seg, skel2area, area_props, len_k, br_num_k, br_pts_k, fwhm_k, dm_k, soam_rad_k, branch_len_k);
        L23_S = aggregate_for_selection(idx_l23_small, seg_col_range{2}, A, top_view_seg, skel2area, area_props, len_k, br_num_k, br_pts_k, fwhm_k, dm_k, soam_rad_k, branch_len_k);

        %% ===== Write to Summary Table =====
        new_row = { ...
            current_tertiary, ...
            vessel_density_global, mean_intensity_global, vessel_area_total, vessel_count, mean_width_global, total_branch_num, total_branch_points, mean_curve_dm_global, mean_soam_rad_global, sum_curve_dm_global, sum_soam_rad_global, max_len_global, mean_branch_len_global, ...
            F13_L.Density,F13_L.Intensity,F13_L.Area,F13_L.Count,F13_L.WidthFWHM,F13_L.BranchNum,F13_L.BranchPoints,F13_L.CurveDM,F13_L.SOAM_Rad,F13_L.sumDM,F13_L.sumSOAM,F13_L.MaxLen,F13_L.MeanBranchLen, ...
            F13_S.Density,F13_S.Intensity,F13_S.Area,F13_S.Count,F13_S.WidthFWHM,F13_S.BranchNum,F13_S.BranchPoints,F13_S.CurveDM,F13_S.SOAM_Rad,F13_S.sumDM,F13_S.sumSOAM,F13_S.MaxLen,F13_S.MeanBranchLen, ...
            L23_L.Density,L23_L.Intensity,L23_L.Area,L23_L.Count,L23_L.WidthFWHM,L23_L.BranchNum,L23_L.BranchPoints,L23_L.CurveDM,L23_L.SOAM_Rad,L23_L.sumDM,L23_L.sumSOAM,L23_L.MaxLen,L23_L.MeanBranchLen, ...
            L23_S.Density,L23_S.Intensity,L23_S.Area,L23_S.Count,L23_S.WidthFWHM,L23_S.BranchNum,L23_S.BranchPoints,L23_S.CurveDM,L23_S.SOAM_Rad,L23_S.sumDM,L23_S.sumSOAM,L23_S.MaxLen,L23_S.MeanBranchLen ...
            };
        results_summary = [results_summary; new_row];

        clear top_view_seg A ve_skel ve_skel_seg L_skel L_area ve_props area_props ...
            branch_points_map len_k area_k br_pts_k br_num_k fwhm_k dm_k soam_rad_k ...
            centroid_cols skel2area
    end
    fprintf('\nProcessed primary folder %d/%d! Time: %.2f seconds\n', primary_idx, total_primary, toc);
end

%% Save Excel Summary
if ~isempty(results_summary)
    excel_file = fullfile(root_path, 'vessel_analysis_summary_3D_25_12_3.xlsx');
    writetable(results_summary, excel_file);
    fprintf('\nResults summary table saved to: %s\n', excel_file);
else
    warning('No result data generated, unable to create Excel summary table');
end