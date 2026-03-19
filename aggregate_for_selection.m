function result = aggregate_for_selection(idx_select, seg_range, A, top_view_seg, skel2area, area_props, ...
    len_k, br_num_k, br_pts_k, fwhm_k, dm_k, soam_rad_k, branch_len_k)

    if isempty(idx_select)
        % If no vessels are found, return NaN or 0 directly
        result = struct('Density', 0, 'Intensity', 0, 'Area', 0, 'Count', 0, ...
                        'WidthFWHM', NaN, 'BranchNum', 0, 'BranchPoints', 0, ...
                        'CurveDM', NaN, 'SOAM_Rad', NaN, 'sumDM', NaN, 'sumSOAM', NaN, ...
                        'MaxLen', NaN, 'MeanBranchLen', NaN);
        return;
    end

    % Vessel count
    count_vessel = numel(idx_select);

    % Area
    area_vals = zeros(1,count_vessel);
    for k = 1:count_vessel
        if skel2area(idx_select(k))>0
            area_vals(k) = numel(area_props(skel2area(idx_select(k))).PixelIdxList);
        else
            area_vals(k) = 0;
        end
    end
    total_area = sum(area_vals);

    % Vessel density (relative to the area of the segmented region)
    seg_mask = false(size(A));
    seg_mask(:, seg_range) = true;
    seg_area_total = nnz(seg_mask);
    density_val = total_area / seg_area_total;

    % Mean vessel intensity (average)
    intensity_vals = zeros(1,count_vessel);
    for k = 1:count_vessel
        if skel2area(idx_select(k))>0
            pix_idx = area_props(skel2area(idx_select(k))).PixelIdxList;
            intensity_vals(k) = mean(double(top_view_seg(pix_idx)));
        else
            intensity_vals(k) = 0;
        end
    end
    mean_intensity = mean(intensity_vals,'omitnan');

    % FWHM, DM, SOAM
    width_vals = fwhm_k(idx_select);
    dm_vals    = dm_k(idx_select);
    soam_vals  = soam_rad_k(idx_select);

    % Branches
    branch_num_vals  = br_num_k(idx_select);
    branch_pts_vals  = br_pts_k(idx_select);

    % Maximum length
    max_len_val = max(len_k(idx_select));

    % Mean branch length
    mean_branch_len_val = mean(branch_len_k(idx_select), 'omitnan');

    % Return results
    result = struct();
    result.Density     = density_val;
    result.Intensity   = mean_intensity;
    result.Area        = total_area;
    result.Count       = count_vessel;
    result.WidthFWHM   = mean(width_vals,'omitnan');
    result.BranchNum   = sum(branch_num_vals);
    result.BranchPoints= sum(branch_pts_vals);
    result.CurveDM     = mean(dm_vals,'omitnan');
    result.SOAM_Rad    = mean(soam_vals,'omitnan');
    result.sumDM       = sum(dm_vals,'omitnan');
    result.sumSOAM     = sum(soam_vals,'omitnan');
    result.MaxLen      = max_len_val;
    result.MeanBranchLen = mean_branch_len_val;
end