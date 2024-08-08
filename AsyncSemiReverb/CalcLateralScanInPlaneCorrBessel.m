function bessel = CalcLateralScanInPlaneCorrBessel(x,k)
  
  % Computes semi-reverberant elastography in-plane Bessel function given a
  % wavenumber k and x array. 
  %
  % 
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute 
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  bessel =  3 * (...
    sin(k .* x) ./ (k .* x) -...
    sin(k .* x) ./ (k .* x) .^ 3 * 2 +...
    cos(k .* x) ./ (k .* x) .^ 2 * 2);
end

