function [k, mseMin, acFit] = KFitLateralScanningInPlaneExcitation...
     (xCorrAxis, displacementXcorr, kMinSearch, kMaxSearch, numSearch)

  % BESSEL K FITTING FOR IN PLANE EXCITATION
  % COMPUTES MEAN SQUARE ERROR FOR ALL K WITHIN THE GIVEN RANGE AND PICKS MIN
  % xCorrAxis:         with dimensions [1, nX] must be 2D 
  % displacementXcorr: with same dimensions as xCorrAxis
  % kMinSearch:        smalled k to consider
  % kMaxSearch:        largest k to consider
  % numSearch:         number of k values to consider, defines search granularity
  % 
  % EXAMPLE
  % [kFit, mse, acFit] = KFitLateralScanningInPlaneExcitation(xAxisFitting, ...
  %     real(acDeltaXOneSided), 0, 15, 1);
  %
  % 
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute 
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
    
  BesselX = xCorrAxis;
  BesselY = displacementXcorr;
  kVec = linspace(kMinSearch, kMaxSearch, numSearch);
  kVec(kVec == 0) = []; % Omit 0 because we can't divid by zero.
  mseVec = CalcLateralScanInPlaneCorrBesselMSE(kVec(:), BesselX(:).', BesselY(:).');
  [mseMin, kIdx] = min(mseVec, [], 1, 'omitnan');
  % Check if this index is the first or last: that means it's not a real minimum
  if kIdx == 1 || kIdx == numel(kVec)
    k = nan;
    mseMin = nan;
  else
    k = kVec(kIdx);
  end
  acFit = CalcLateralScanInPlaneCorrBessel(BesselX, k);
  
end

