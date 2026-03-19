function vesselStats = analyzeVesselsMultiFrame(mask1, fig1)
% analyzeVesselsMultiFrame - Multi-frame vessel mask analysis and segmented statistics
% Input:
%   mask1 - Binary vessel segmentation map [height x width x frames]
%   fig1  - Original image                [height x width x frames]
% Output:
%   vesselStats - Structure containing 36 statistical variables across three segments and total vessel counts

[nRows, nCols, nFrames] = size(mask1);
if size(fig1,1)~=nRows || size(fig1,2)~=nCols || size(fig1,3)~=nFrames
    error('Dimensions of mask1 and fig1 do not match');
end

% Define frame segments
segments = {[3:10],[11:20],[21:30]};

% Initialize output structure
vesselStats = struct();

% Loop through each segment
for s = 1:length(segments)
    segFrames = segments{s};
    
    % Initialize statistical arrays for each frame
    small_num = zeros(length(segFrames),1);
    medium_num = zeros(length(segFrames),1);
    large_num = zeros(length(segFrames),1);
    
    small_area = zeros(length(segFrames),1);
    medium_area = zeros(length(segFrames),1);
    large_area = zeros(length(segFrames),1);
    
    small_strength = zeros(length(segFrames),1);
    medium_strength = zeros(length(segFrames),1);
    large_strength = zeros(length(segFrames),1);
    
    small_width = zeros(length(segFrames),1);
    medium_width = zeros(length(segFrames),1);
    large_width = zeros(length(segFrames),1);
    
    % Process each frame
    for f = 1:length(segFrames)
        frameIdx = segFrames(f);
        maskFrame = mask1(:,:,frameIdx);
        imgFrame = fig1(:,:,frameIdx);
        
        % --- Core Correction: Use bwconncomp to separate disconnected vessels ---
        CC = bwconncomp(maskFrame);
        stats = regionprops(CC, 'Area', 'EquivDiameter', 'MajorAxisLength', ...
            'MinorAxisLength', 'PixelIdxList');
        
        % Classify vessels
        smallVessels = [];
        mediumVessels = [];
        largeVessels = [];
        smallAreas = [];
        mediumAreas = [];
        largeAreas = [];
        
        for i = 1:length(stats)
            area = stats(i).Area;
            if area >= 4 && area <= 49
                smallVessels = [smallVessels; i];
                smallAreas = [smallAreas; area];
            elseif area >= 50 && area <= 99
                mediumVessels = [mediumVessels; i];
                mediumAreas = [mediumAreas; area];
            elseif area > 99
                largeVessels = [largeVessels; i];
                largeAreas = [largeAreas; area];
            end
        end
        
        % Calculate binary map and intensity for each vessel class
        classes = {smallVessels, mediumVessels, largeVessels};
        areasList = {smallAreas, mediumAreas, largeAreas};
        widthsList = {[],[],[]};
        strengths = zeros(1,3);
        
        for cls = 1:3
            indices = classes{cls};
            classAreas = areasList{cls};
            bwClass = false(size(maskFrame));
            for idx = indices'
                bwClass(stats(idx).PixelIdxList) = 1;
            end
            
            % Vessel intensity
            multiplied = abs(double(imgFrame).*bwClass);
            if any(multiplied(:)>0)
                strengths(cls) = mean(multiplied(multiplied>0));
            else
                strengths(cls) = 0;
            end
            
            % Vessel width
            widths = zeros(length(indices),1);
            for j = 1:length(indices)
                AR = stats(indices(j)).MinorAxisLength / stats(indices(j)).MajorAxisLength;
                if AR >= 0.5
                    widths(j) = stats(indices(j)).EquivDiameter / 2;
                else
                    vesselMask = false(size(maskFrame));
                    vesselMask(stats(indices(j)).PixelIdxList) = 1;
                    skel = bwmorph(vesselMask,'skel',Inf);
                    skel = bwmorph(skel,'spur',5); % Remove short spurs
                    skelLength = nnz(skel);
                    if skelLength>0
                        widths(j) = stats(indices(j)).Area / skelLength;
                    else
                        widths(j) = 0;
                    end
                end
            end
            widthsList{cls} = widths;
        end
        
        % Save per-frame data
        small_num(f) = length(smallVessels);
        medium_num(f) = length(mediumVessels);
        large_num(f) = length(largeVessels);
        
        % Area
        if ~isempty(smallAreas)
            small_area(f) = sum(smallAreas,'omitnan');
        else
            small_area(f) = 0;
        end
        if ~isempty(mediumAreas)
            medium_area(f) = sum(mediumAreas,'omitnan');
        else
            medium_area(f) = 0;
        end
        if ~isempty(largeAreas)
            large_area(f) = sum(largeAreas,'omitnan');
        else
            large_area(f) = 0;
        end

        small_strength(f) = strengths(1);
        medium_strength(f) = strengths(2);
        large_strength(f) = strengths(3);
        
        % Width
        if ~isempty(widthsList{1})
            small_width(f) = mean(widthsList{1},'omitnan');
        else
            small_width(f) = 0;
        end
        if ~isempty(widthsList{2})
            medium_width(f) = mean(widthsList{2},'omitnan');
        else
            medium_width(f) = 0;
        end
        if ~isempty(widthsList{3})
            large_width(f) = mean(widthsList{3},'omitnan');
        else
            large_width(f) = 0;
        end
    end
    
    % Average over segments and store results
    vesselStats(s).segment = s;
    vesselStats(s).small_num = mean(small_num);
    vesselStats(s).medium_num = mean(medium_num);
    vesselStats(s).large_num = mean(large_num);
    
    vesselStats(s).small_area = mean(small_area);
    vesselStats(s).medium_area = mean(medium_area);
    vesselStats(s).large_area = mean(large_area);
    
    vesselStats(s).small_strength = mean(small_strength);
    vesselStats(s).medium_strength = mean(medium_strength);
    vesselStats(s).large_strength = mean(large_strength);
    
    vesselStats(s).small_width = mean(small_width);
    vesselStats(s).medium_width = mean(medium_width);
    vesselStats(s).large_width = mean(large_width);
    
    % Calculate total vessel count for the segment
    vesselStats(s).total_num = sum(small_num + medium_num + large_num);
end

end