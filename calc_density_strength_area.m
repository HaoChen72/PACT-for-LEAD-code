function [density_seg1, density_seg2, density_seg3, ...
          intensity_seg1, intensity_seg2, intensity_seg3, ...
          area_seg1, area_seg2, area_seg3] = calc_density_strength_area(resized_masks, resized_labels, resized_figs)
% Calculate vessel density, signal intensity, and vessel area for three segments (3-10, 11-20, 21-30)
% Output: 9 independent variables

    % Define the range for the three segments
    segments = {3:10, 11:20, 21:30};
    nSeg = numel(segments);

    % Temporary storage for results
    density_mean = zeros(1, nSeg);
    intensity_mean = zeros(1, nSeg);
    area_mean = zeros(1, nSeg);

    for s = 1:nSeg
        frames = segments{s};
        numFrames = numel(frames);

        density_vals = zeros(1, numFrames);
        intensity_vals = zeros(1, numFrames);
        area_vals = zeros(1, numFrames);

        for i = 1:numFrames
            f = frames(i);

            % mask, label, fig
            mask_f = resized_masks(:, :, f);
            label_f = resized_labels(:, :, f);
            fig_f = resized_figs(:, :, f);

            % Combined mask = current frame's mask
            combined_mask = (mask_f == 1);

            % Vessel area
            area_vals(i) = nnz(combined_mask);

            % Tissue area (total foot pixels)
            foot_pixel = nnz(label_f);

            % Vessel density
            if foot_pixel == 0
                density_vals(i) = NaN;
            else
                density_vals(i) = area_vals(i) / foot_pixel;
            end

            % Blood signal intensity
            if area_vals(i) == 0
                intensity_vals(i) = NaN;
            else
                intensity_vals(i) = mean(mean(abs(fig_f(combined_mask))));
            end
        end

        % Mean value for each segment
        density_mean(s)   = mean(density_vals, 'omitnan');
        intensity_mean(s) = mean(intensity_vals, 'omitnan');
        area_mean(s)      = mean(area_vals, 'omitnan');
    end

    % Split into 9 independent variables
    density_seg1   = density_mean(1);
    density_seg2   = density_mean(2);
    density_seg3   = density_mean(3);

    intensity_seg1 = intensity_mean(1);
    intensity_seg2 = intensity_mean(2);
    intensity_seg3 = intensity_mean(3);

    area_seg1      = area_mean(1);
    area_seg2      = area_mean(2);
    area_seg3      = area_mean(3);
end