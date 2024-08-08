function [phiUnwrapped] = Perform2DDCTPhaseUnwrapping(psi)
  % Perform2DDCTPhaseUnwrapping 2D phase unwrapping following:
  %    Ghiglia DC, Romero LA. Robust two-dimensional weighted and unweighted
  %    phase unwrapping that uses fast transforms and iterative methods. J Opt
  %    Soc Am A, JOSAA. Optica Publishing Group; 1994 Jan 1;11(1):107–117.
  % Also inspired by Firman (2024). 2D Weighted Phase Unwrapping 
  %    (https://www.mathworks.com/matlabcentral/fileexchange/60345-2d-weighted-phase-unwrapping), 
  %    MATLAB Central File Exchange. 
  % Uses the discrete cosine transform to find the LSE solution in 2D
  % 
  %
  % 
  % Authors:  Ginger Schmidt (1,2), Néstor Uribe-Patarroyo (1) 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA
  % 2. Institute for Medical Engineering and Science, Massachusetts Institute 
  % of Technology, 77 Massachusetts Avenue, Cambridge,, MA, USA
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  
  % Compute the wrapped wrapped-phase differences along 1st index
  deltaZ = wrapToPi(diff(psi, 1, 1));
  % Set boundary conditions
  deltaZ = padarray(deltaZ, [1 0], 'both');
  % Compute the wrapped wrapped-phase differences along 2nd index
  deltaX = wrapToPi(diff(psi, 1, 2));
  % Set boundary conditions
  deltaX = padarray(deltaX, [0 1], 'both');
  
  % Calculate rho, essentially the known side of the Poisson equation (the
  % second derivatives of psi)
  rho = diff(deltaZ, 1, 1) + diff(deltaX, 1, 2);

  % Get the unwrapped phase by solving the Poisson equation
  phiUnwrapped = SolvePoissonWithDCT(rho);
end

function phi = SolvePoissonWithDCT(rho)
  % Solve the Poisson equation using dct
  % Apply 2D DCT
  rhoDCTZX = dct(dct(rho, [], 1), [], 2);
  % Calculate vectors and factors needed for the solution. We want these to be
  % in double precision
  [nZ, nX] = size(rho, [1 2]);
  [xMat, zMat] = meshgrid(0:nX - 1, 0:nZ - 1);
  % Now the solution is the 2D DCT multiplied by the corresponding factors
  phiDCTZX = rhoDCTZX ./ (2 .* (cos(pi .* xMat ./ nX) + cos(pi .* zMat ./ nZ) - 2));
  % We preserve the bias (mean) of the wrapped phase by making the undefined 0,0
  % term equal to rhoDCTZX(0,0) (instead of 0 bias like the old function).
  phiDCTZX(1, 1, :) = rhoDCTZX(1, 1, :);
  % Now apply inverse DCT to get the result
  phi = idct(idct(phiDCTZX, [], 1), [], 2);
end