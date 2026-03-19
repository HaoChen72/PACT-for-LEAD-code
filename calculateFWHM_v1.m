function [FWHM, data] = calculateFWHM_v1(data, varargin)
  if nargin > 1
    gaussParam = varargin{1};
  else
    gaussParam = 0;
  end
  if length(data) >= 9
     data = interp(data, 2);
  end
  if gaussParam
    data = imgaussfilt(data, gaussParam);
  end
  FWHM = NaN;
  
  [peakVal, peakIdx] = findpeaks(data, 'MinPeakHeight', max(data)*0.5, 'WidthReference', 'halfheight');
  [peakVal, peakIdxVal] = max(peakVal);
  peakIdx = peakIdx(peakIdxVal);
  
  
  if isempty(peakIdx)
    peakIdx = NaN;
    return;
  end
  
  halfMax = peakVal(1) *0.5;
  leftIdx = find(data(1:peakIdx(1)) <= halfMax, 1, 'last');
  rightIdx = find(data(peakIdx(1):end) <= halfMax, 1, 'first') + peakIdx(1) - 1;
  
  if isempty(leftIdx) || leftIdx == 1
    leftIdx = NaN;
  end
  if isempty(rightIdx) || rightIdx == length(data)
    rightIdx = NaN;
  end
  if ~isnan(leftIdx) && ~isnan(rightIdx)
    FWHM = rightIdx - leftIdx;
  end

  
end

function filteredArray = GaussianFilter(array, sigma)
  filterSize = 2 * ceil(3 * sigma) + 1;
  gaussFilter = exp(-((1:filterSize) - ceil(filterSize / 2)).^2 / (2 * sigma^2));
  gaussFilter = gaussFilter / sum(gaussFilter);
  filteredArray = conv(array, gaussFilter, 'same');
end
