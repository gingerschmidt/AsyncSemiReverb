% To process your own data, you should provide a variable tom (supplied
% here by tom.mat) that contains a complex-valued tomogram with dimensions [Z, X, Bscans,
% Polarization Channels] AND update other parameters.
load tom.mat
tom = tomSaved;

%% ---------------------------------------------------------------------
% Define data parameters
% MODIFY THIS SECTION FOR YOUR OWN DATA
% ----------------------------------------------------------------------
% DEFINE tomogram parameters
nAlinesPerBscanNominal =...
  1792;                  % *BEFORE* any cropping, from original data (important for demodulation)
nBscanStepSize = 2;      % Need at least 2 B-scans at each Y location to measure displacement
scanWidthX = 14;         % in mm; fast-axis galvo scan range
scanWidthY = 14;         % in mm; slow-axis galvo scan range
noiseFloorROI =...
  {24:32, 500:700};      % Define a region {z0:zEnd, x0:xEnd} in air to calculate the noise floor for intensity masking
xROI = 277:1792;         % For cropping tomogram, omitting galvo flyback
zROI = 1:900;            % For cropping tomogram, omitting depth areas with no signal

% Displacement calculation parameters
dispFilterRadii =...
  [4, 0];                % in px; size of Gaussian kernel radii for spatially filtering the Doppler OCT displacement

% Elastography parameters
excitationFreq = 1981;   % in Hz; chosen to maximize displacement between B-scans
alineRate = 100000;      % in Hz; OCT laser A-line repetition rate
filtWidth = 30;          % in px; Fourier-domain filter size, wide enough to capture displacement signal but avoid overlapping with DC.
refractiveIdx = 1.4;     % Refractive index in tissue (approx.)  for Doppler OCT displacement calculation
wavelengthOCT = 1.3;     % in µm for Doppler OCT displacement calculation
corrWindowHeight = 10;   % in px; Correlation window size in depth
dimAC = 2;               % Autocorrelation dimension (2 is x, code only works with this value)
aveDimsAC = [1 3];       % Other dimensions averaged for the autocorrelation calculation (1 is z and 3 is y)
nZeroCrossing = 3;       % Number of zero-crossing used for fitting k
numSearch = 1001;        % Number k values to search, defines search granularity
kMinSearch = 0;          % in µm^{-1}; min k value to consider
kMaxSearch = 15;         % in µm^{-1}; max k value to consider
    
% Visualization & graphing parameters 
currFig = 0;             % Code generates figure numbers offset by currFig. Useful to compare different runs by changing currFig to 100, 200,...
logLim = [70, 120];      % in dB; limits for tomogram intensity plotting, adjust for your tomogram dynamic range
displacementLim =...
  [-0.5, 0.5];           % in µm; for displacement plotting
nYToPlot = 10;           % Number of y-locations to visualize (< nY)
nDispToVis = 10;         % Number of displacement frames to visualize (< nBscans/2)
pauseTime = 0.1;         % in s; rate of refresh for visualizations
maskThreshold = 9.5;     % in dB; tune this if the YZ view intensity mask isn't right.

% Surface detection and flattening parameters
surfaceInterfaceDb = 89; % in dB; tune this value if surface detection is failing
bigRegionFactor = 7;     % >= 3 for uniformly strong surfaces, lower for surfaces with breaks
maxSurfaceGradient = 15; % in px/px; used to detect artifactual abrupt changes in tissue height
dcMask = [];             % if your tomogram has a strong DC artifact, provide range so it is removed during surface detection
filterKernel =...
  ones([20, 30, 1]);     % in px; kernel used for filtering the detected surface. Tune this value if surface detection is failing

% Generate displacement at *any* time step from just 2 Bscans
nPeriods = 2;            % Number of periods to generate
nTimeStepsPerPeriod =...
  41;                    % Number of time steps to generate, per period
generateThisY = 52;      % Generate displacements for this y location (< nY)

%% ---------------------------------------------------------------------
% Compute other data parameters
% ----------------------------------------------------------------------
% Add subfolders to path 
addpath(genpath('SupportFunctions'));
addpath(genpath('AsyncSemiReverb'));

% Crop tomogram
tom = tom(zROI, xROI, :, :);
% Calculate noise floor
noiseFloorDb = 10 * log10(mean(abs(tom(noiseFloorROI{:}, :, :)) .^ 2, [1 2 3]));

% Shear wave excitation and demodulation parameters.
% The demodulation carrier frequency is a function excitation frequency, number
% of samples, and the laser sweep repition rate. [See Section 2.1.2]
[nZ, nX, nBscans, nPolChannels] = size(tom);
nY = nBscans / nBscanStepSize;
demodulationShiftPx = excitationFreq * nX / alineRate;
% For accurate demodulation, we need this value to be very close to an
% integer pixel
fprintf('Fractional part of demodulationShiftPx: %g\n', demodulationShiftPx - round(demodulationShiftPx))
if abs(demodulationShiftPx - round(demodulationShiftPx)) > 0.1
  warning('Are you *sure* you want to proceed with such large fractional demodulation px?')
  % To get a more round demodulation pixel, try cropping the tomogram along X to
  % change the FFT length.
end
demodulationShiftPx = round(demodulationShiftPx);
vScan = scanWidthX / nAlinesPerBscanNominal * alineRate;  % in mm/sec
fs = nAlinesPerBscanNominal / scanWidthX;                 % in 1/mm
scanningFreq = excitationFreq / vScan;                    % in 1/mm
dopplerFactor = wavelengthOCT / (4 * pi * refractiveIdx); % in µm/rad; to convert from phase difference to µm

% Visualization & graphing Parameters 
LATEX_DEF = {'Interpreter', 'latex'};
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
diffCMap = colorcet('C3');
cMapDisp = colorcet('D6');
nanColor = 0.5 * [1, 1, 1];

%% ---------------------------------------------------------------------
% Surface detection & flattening
% ----------------------------------------------------------------------
% It is important to flatten the tissue since we are computing the
% autocorrelation across the entire x dimension. For layered samples, this
% ensures that we do not average across different tissue types. It also ensures
% we are not averaging tissue with air. Flattening may not always be neccessary.

tissueFlatteningOptions = struct('surfaceInterfaceDb', surfaceInterfaceDb,...
  'bigRegionFactor', bigRegionFactor, 'dcMask', dcMask,...
  'nBscanStepSize', nBscanStepSize, 'filterKernel', filterKernel,...
  'desiredTissueSurfaceDepth', [], 'maxSurfaceGradient', maxSurfaceGradient,...
  'figNum', currFig + 1);

[tomFlatX, targetSurfaceZPerY, meanSurfaceZPerY] = DetectTissueSurfaceAndFlatten(tom, tissueFlatteningOptions);

%% ---------------------------------------------------------------------
% Compute displacement, unwrap, & surface wave correction
% ----------------------------------------------------------------------
% Compute displacement from tomogram (in um) with Doppler OCT. 
% We assume tissue all has the same refractive index.
[displacementRaw, displacementPhasorRaw] = ...
  CalcDisplacementFromTom(tomFlatX, dopplerFactor,...
  'bscanGroupSize', nBscanStepSize, 'doSpatialFiltering', true, ...
  'filterRadii', dispFilterRadii);
% Unwrap phase, following Ghiglia DC and Romero LA, Optica Publishing Group (1994)
fprintf('Unwrapping %d Bscans...', size(displacementRaw, 3));
displacementUnwrappedX = (Perform2DDCTPhaseUnwrapping(displacementRaw / dopplerFactor) * dopplerFactor);
fprintf('done.\n');

% Surface wave correction, following Song et. al., J. Biomed. Opt 18, 121505 (2013). 
fprintf('Applying surface wave correction...');
displacementBscanStepSize = nBscanStepSize * (nBscanStepSize - 1) / nBscanStepSize;
displacementSurfWaveCorr = displacementUnwrappedX;
for thisY = 1:nY
  surfaceWaveOffset = (refractiveIdx - 1) *...
    displacementUnwrappedX(targetSurfaceZPerY(:, thisY), :, displacementBscanStepSize * (thisY - 1) + (1:displacementBscanStepSize));
  displacementSurfWaveCorr(:, :, displacementBscanStepSize * (thisY - 1) + (1:displacementBscanStepSize)) =...
    displacementUnwrappedX(:, :, displacementBscanStepSize * (thisY - 1) + (1:displacementBscanStepSize)) + surfaceWaveOffset;
end
fprintf('done.\n');
% Unwrap phase again after surface wave correction
fprintf('Unwrapping %d Bscans after surface wave correction...', size(displacementRaw, 3));
displacementSurfWaveCorr = (Perform2DDCTPhaseUnwrapping(displacementSurfWaveCorr / dopplerFactor) * dopplerFactor);
fprintf('done.\n');

%% ---------------------------------------------------------------------
% Demodulation for coherent wave recovery
% ----------------------------------------------------------------------
% In this section, we demodulate the effect of raster-scanning to recover the
% coherent displacement field. 
% [Refer to sections 2.1.2 and 2.3 for further details.]
% First, compute the Fourier transform of the displacement and average 
% the power spectrum across time steps. 
displacementScanFTX = displacementSurfWaveCorr;
displacementScanFTX =  fftshift(fft(fftshift(displacementScanFTX, 2), [], 2), 2);
displacementScanPS = mean(abs(displacementScanFTX) .^ 2, 3);
% Then, we downshift (with circshift) and bandpass filter. 
% Finally, compute the inverse Fourier transform.
displacementDownshiftedFTX = displacementScanFTX;
displacementDownshiftedFTX = circshift(displacementDownshiftedFTX, demodulationShiftPx, 2);
displacementDownshiftedFTX(:, [1:end / 2 + 1 - filtWidth / 2, end / 2 + 1 + filtWidth / 2:end], :) = 0;
displacementDemod = ifftshift(ifft(ifftshift(displacementDownshiftedFTX, 2), [], 2), 2);

%% ---------------------------------------------------------------------
% Generate displacement at *any* time step from just 2 Bscans
% ----------------------------------------------------------------------
% As shown in Supplemental Video 1, we can generate the coherent displacement
% field at any time step from just 2 B-scans [see Section 2.1.2].
% Multiply the complex displacement by a corresponding complex phase to generate
% the displacement at a different point in time. 
% Compute the real part of the displacement for visualization. 
displacementPerYStep = reshape(displacementDemod, nZ, nX, displacementBscanStepSize, nY);
thisPhaseVec = linspace(0, 2 * pi * nPeriods, nPeriods * nTimeStepsPerPeriod);
for thisPhase = thisPhaseVec
  figure(currFig + 2)
  imagesc(squeeze(real(...
    displacementPerYStep(max(targetSurfaceZPerY(generateThisY) - 50, 1):end, :, :, generateThisY)...
    .* exp(-1i .* thisPhase))), displacementLim), colorbarlbl('[\mum]'), colormap(cMapDisp)
  % Modify displacementLim above if plotting is out of range!
  title('Demodulated displacement at arbitrarty time step'), xlabel('x'), ylabel('z')
  drawnow
end

%% ---------------------------------------------------------------------
% Asynchronous Semi-Reverberant Elastography
% ----------------------------------------------------------------------
% At each y-location, we compute the 1D autocorrelation across the entire x
% dimension. The code may be adapted for different correlation window sizes by
% changing corrWindowSizeDeltaX.
% At each y-location, we compute k across depth by fitting a Bessel function sum
% to the autocorrelation. We use the autocorrelation only up to a specified zero
% crossing to avoid overfitting in some regions. 
% [Refer to Section 2.1.3 and Equation 4 for further details.] 
corrWindowSizeDeltaX = [corrWindowHeight, width(displacementUnwrappedX)];
if mod(corrWindowSizeDeltaX,2) == 1
  warning('Please redfine your tomogram so the correlation dimension (nX) is even.')
end
nWindowsZ = floor(height(displacementUnwrappedX) / corrWindowSizeDeltaX(1));

% Determine k Bessel fitting parameters
meanSubtractionMethod = 'ACWindow'; % Mean subtraction can be 'none' or 'ACWindow'
dk = 1 / (scanWidthX * 1E3);        % in µm^-1
kAxis = dk * (-corrWindowSizeDeltaX(2):corrWindowSizeDeltaX(2));
xAxisSpacing = scanWidthX / nX;     % in mm/px
kFitResult = zeros(nWindowsZ, nY);  % Allocate empty arrays
kFitResultSSE = zeros(nWindowsZ, nY); % Allocate empty arrays
fprintf('Processing Y position (out of %d): ', nY)
tic
for thisY = 1:nY
  fprintf('%d, ', thisY)
  displacementThisY = displacementDemod(:, :, displacementBscanStepSize * (thisY - 1) + 1);
  % Calculate 1D autocorrelations
  [acDeltaX, displacementACFT] = CalcWindowed1DAutocorrelation(...
    displacementThisY, corrWindowSizeDeltaX, dimAC, aveDimsAC, meanSubtractionMethod);
  deltaXZero = corrWindowSizeDeltaX(2);
  % Normalize ACs
  acDeltaX = squeeze(acDeltaX).';
  acDeltaX = acDeltaX ./ abs(acDeltaX(:, deltaXZero));
  acDeltaXOneSided = acDeltaX(:, deltaXZero:end);
  % K fitting
  kFitThisY = zeros(nWindowsZ, 1, 'single');
  mseThisY = zeros(nWindowsZ, 1, 'single');
  acFitThisY = zeros(nWindowsZ, corrWindowSizeDeltaX(2), 'single');
  for thisZWin = 1:nWindowsZ
    % Zero crossing condition is a pixel in which the array changes sign
    negativeACIdx = find(real(acDeltaXOneSided(thisZWin, 2:end)) .* real(acDeltaXOneSided(thisZWin, 1:end - 1)) < 0,...
      nZeroCrossing, 'first');
    if numel(negativeACIdx) == 0
      % No zero crossing found
      fitIdx = round(3 * corrWindowSizeDeltaX(2) / 4);
    elseif numel(negativeACIdx) >= nZeroCrossing
      % Desired zero crossing found
      negativeACIdx = negativeACIdx(nZeroCrossing);
      fitIdx = min(negativeACIdx, corrWindowSizeDeltaX(2) - 1);
    else
      % Other zero crossings found, pick closest to our desired one
      negativeACIdx = negativeACIdx(end);
      fitIdx = min(negativeACIdx, corrWindowSizeDeltaX(2) - 1);
    end
    xAxisFitting = (1:fitIdx) * xAxisSpacing;
    [kFitThisY(thisZWin), mseThisY(thisZWin), acFitThisY(thisZWin, 2:fitIdx + 1)] = ...
      KFitLateralScanningInPlaneExcitation(xAxisFitting, ...
      real(acDeltaXOneSided(thisZWin, 2:fitIdx + 1)), ...
      kMinSearch, kMaxSearch, numSearch);
  end
  % Save result
  kFitResult(:, thisY) = kFitThisY;
  kFitResultSSE(:, thisY) = mseThisY;
end
elapsedTime = toc;
fprintf('done (took %.1f seconds).\n', elapsedTime)

%% ---------------------------------------------------------------------
% Asynchronous Semi-Reverberant Elastography RESULTS
% ----------------------------------------------------------------------

% Format tomogram for YZ viewing and compute intensity mask
% If the mask doesn't look right, tune maskThreshold.
tomACMeanInt = squeeze(mean(reshape(mean(sum(abs(tomFlatX) .^ 2, 4), dimAC),...
  corrWindowSizeDeltaX(1), nWindowsZ, nBscanStepSize, nY), aveDimsAC));
tomACMeanIntMask = 10 * log10(tomACMeanInt) > maskThreshold + mean(noiseFloorDb);
for thisY = 1:nY
  tomACMeanIntMask(1:floor(targetSurfaceZPerY(thisY) / corrWindowSizeDeltaX(1)), thisY) = false;
end
tomACMeanIntMask = single(imclose(tomACMeanIntMask, strel('disk', 5, 0)));
tomACMeanIntMask(~tomACMeanIntMask) = nan;

% Calculate shear modulus from k [mm^-1]
shearSpeedMapFit = 2 * pi * excitationFreq ./ kFitResult; % mm/s
materialDensity = 1E-3; % g/mm^3 (1E3 kg/m^3=1E3 kg/mm^3/1E9=1E3µg/mm^3), assuming same as water
shearModulusFit = shearSpeedMapFit .^ 2 * materialDensity; % Pa (mm^2/s^2*g/mm^3=g/(mm*s^2)=1E-3kg/(1E-3m*s^2)=kg/(ms^2)=Pa)

%% Plot results
figure(currFig + 3)
imagesc(10 * log10(tomACMeanInt) .* tomACMeanIntMask, logLim)
title('\bf Tomogram Intensity (YZ view)'), colorbarlbl('[dB]'), colormap(gray(256))
xlabel('y'), ylabel('z')
figure(currFig + 4)
subplot(2,1,1)
imagescnan(tomACMeanIntMask .* kFitResult, [1 5])
title('\bf Wavenumber (k) fit in [1/mm]'), colorbarlbl('[1/mm]'), colormap(viridis)
xlabel('y'), ylabel('z')
subplot(2,1,2)
imagescnan(tomACMeanIntMask .* shearModulusFit, [3E3 9E4])
title('\bf Shear modulus (G) fit in [Pa]'), colorbarlbl('[Pa]'), colormap(viridis)
xlabel('y'), ylabel('z')
hAx = gca; hAx.ColorScale = 'log';



