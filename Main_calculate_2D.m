%% Automated Vessel Parameter Calculation and Excel Export Script
% Author: Dong Shirui
% Date: 2025-08-22
% Input: Set root_path

%% Initialization Settings
root_path = 'K:\footData\'; % Modify to your data path
clc;
close all force;
warning('off','images:initSize:adjustingMag'); % Disable image-related warnings

%% Path Validation and Directory Tree Construction
if ~isfolder(root_path)
    error('Root path does not exist: %s', root_path);
end

% Get list of primary subfolders (Patient IDs)
primary_list = dir(fullfile(root_path, '*'));
primary_list = primary_list([primary_list.isdir]);
primary_list = primary_list(~ismember({primary_list.name}, {'.','..'}));
total_primary = numel(primary_list);
fprintf('Found %d primary folders\n', total_primary);

%% Create Results Summary Table
results_summary = table();

%% Main Processing Flow
for primary_idx = 1:total_primary
    current_primary = fullfile(primary_list(primary_idx).folder,...
                               primary_list(primary_idx).name);
    
    % Locate 3D_str directory in secondary subfolders
    secondary_path = fullfile(current_primary, '3D_str');
    if ~isfolder(secondary_path)
        fprintf('[WARNING] Skipping directory missing 3D_str: %s\n', current_primary);
        continue;
    end
    
    % Get tertiary subfolders L/R (maximum of first two)
    tertiary_list = dir(fullfile(secondary_path, '*'));
    tertiary_list = tertiary_list([tertiary_list.isdir]);
    tertiary_list = tertiary_list(~ismember({tertiary_list.name}, {'.','..'}));
    process_tertiary = tertiary_list(1:min(2,end));
    
    %% Tertiary Folder Processing
    for tertiary_idx = 1:numel(process_tertiary)
        current_tertiary = fullfile(process_tertiary(tertiary_idx).folder,...
                                    process_tertiary(tertiary_idx).name);
        data_file = fullfile(current_tertiary, 'resize_divide.mat');
        if ~isfile(data_file)
            fprintf('[ERROR] Missing data file: %s\n', data_file);
            continue;
        end
        
        %% Load Data
        try
            temp = load(data_file, 'resized_figs','resized_labels','resized_masks');
            figs = temp.resized_figs;
            labels = temp.resized_labels;
            masks = temp.resized_masks;
            
            % Extract shallow regions
            for i = 1:30
                labels1(:,:,i) = logical(erode_mask_n_pix(labels(:,:,i),0) - erode_mask_n_pix(labels(:,:,i),64));
                figs1(:,:,i) = labels1(:,:,i).*figs(:,:,i);
                masks1(:,:,i) = logical(labels1(:,:,i).*masks(:,:,i));
            end
            
            %% Calculate Vessel Density, Intensity, and Area
            [density_seg1, density_seg2, density_seg3,...
             strength_seg1, strength_seg2, strength_seg3,...
             area_seg1, area_seg2, area_seg3] = calc_density_strength_area(masks1, labels1, figs1);

            %% Calculate Large/Medium/Small Vessel Parameters
            vesselStats = analyzeVesselsMultiFrame(masks1, figs1);

            %% Assemble row parameters and add FolderName
            FolderName = current_tertiary;
            rowParams = assembleRowParams(density_seg1,density_seg2,density_seg3,...
                                          strength_seg1,strength_seg2,strength_seg3,...
                                          area_seg1,area_seg2,area_seg3,...
                                          vesselStats, FolderName);

            % Add to summary table
            results_summary = [results_summary; rowParams];

        catch ME
            fprintf('[LOAD/CALCULATION FAILED] %s\nReason: %s\n', data_file, ME.message);
            continue;
        end
    end % End of tertiary folder loop
    fprintf('\nProcessing complete!\nTotal time elapsed: %.2f seconds\n', toc);
end % End of primary folder loop

%% Save Results Summary Table to Excel
if ~isempty(results_summary)
    excel_file = fullfile(root_path, 'vessel_analysis_2D_divide.xlsx');
    writetable(results_summary, excel_file);
    fprintf('\nResults summary table saved to: %s\n', excel_file);
else
    warning('No result data generated, unable to create Excel summary table');
end