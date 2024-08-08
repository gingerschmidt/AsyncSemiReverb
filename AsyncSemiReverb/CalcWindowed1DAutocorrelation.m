function [acArray, arrayPowerSpectrum] = CalcWindowed1DAutocorrelation ...
    (array, windowSize, dim, aveDims, varargin)
  
  % CALCULATE 1D AUTOCORRELATION ALONG GIVEN DIMENSION
  % array:        [Z, X]
  % windowSize:   [windowHeight, X] must be 2D
  % dim:          Dimensions over which we compute the AC
  % aveDims:      Dimensions we average over
  % meanSubtractionMethod: 'none' or 'ACWindow'
  % 
  % EXAMPLE
  % [acDeltaX, displacementACFT] = CalcWindowed1DAutocorrelation(...
  %   displacementThisY, corrWindowSizeDeltaX, dimAC, aveDimsAC, meanSubtractionMethod);
  %
  %
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute 
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  if nargin < 5
    meanSubtractionMethod = 'none';
  else
    meanSubtractionMethod = varargin{1};
  end
  
  % Compute AC dimensions
  arrayDims = size(array);
  nDims = numel(arrayDims);
  colonOp = repmat({':'}, 1, nDims);
  nWindowsZ = floor(arrayDims(1) / windowSize(1));
  nWindowsX = floor(arrayDims(2) / windowSize(2));
  % If not an even multiple, the array will be cropped
  if nWindowsZ * windowSize(1) ~= arrayDims(1) || nWindowsX * windowSize(2) ~= arrayDims(2)
    array = array(1:nWindowsZ * windowSize(1), 1:nWindowsX * windowSize(2), colonOp{3:end});
  end
  
  % Reshape and permute array for convenient processing
  % Window sizes will be in dims 1 and 2, and window indices in dims 3 and 4
  array = reshape(array, windowSize(1), nWindowsZ, windowSize(2), nWindowsX, arrayDims(3:end));
  array = permute(array, [1 3 2 4:nDims + 2]);

  % Subtract mean according to desired method
  if strcmp(meanSubtractionMethod, 'ACWindow')
    array = array - mean(array, [dim aveDims]);
  end
  
  % Autocorrelation
  aveDims(aveDims > 2) = aveDims(aveDims > 2) + 2;
  acDimSize = size(array, dim);
  padSize = zeros(1, nDims);
  padSize(dim) = acDimSize - 1;
  acArray = abs(fft(padarray(array, padSize, 'post'), [], dim)) .^ 2;
  acArray = mean(acArray, aveDims);
  arrayPowerSpectrum = fftshift(acArray, dim);
  acArray = fftshift(ifft(acArray, [], dim), dim);
  
  % Remove bias by dividing by number of actual elements for each displacement
  biasMat = shiftdim(triang(2 * acDimSize - 1), 1 - dim);
  acArray = acArray  ./ biasMat;
end