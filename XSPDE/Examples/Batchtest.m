%   BATCHTEST is a Matlab script to test the operation of an xSPDE toolbox 
%   There are currently 19 projects tested, each with graphical output 
%   Matlab Parallel toolbox is required for a complete test
%   Total runtime should be less than 60s, depending on CPU speed.
%   To see the generated graphs, add '%' before the  'close all' commands

e(1) = Wiener;
e(2) = SHO;
e(3) = Kubo;
e(4) = Soliton;
e(5) = Gaussian;
e(6) = Planar;
e(7) = Gain;
e(8) = Characteristic;
e(9) = Equilibrium;
e(10) = GainTransfer;
close all;         %%Deletes all figures if not wanted
e(11) = GPE;
e(12) = GPE2;
e(13) = Kubotest;
e(14) = SolitonAv;
e(15) = Nonlinear;
e(16) = GaussianAv1D;
e(17) = GaussianAv2D;
e(18) = GaussianAv3D;
e(19) = GaussianKinetic;
close all;         %%Deletes all figures if not wanted
fprintf('xSPDE batch test error score = %f \n',prod(e)*10^45);
fprintf('Expected total error score   = 2.074461 \n' );


