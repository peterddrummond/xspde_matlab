%   BATCHTEST2 is a Matlab script to test the operation of the xSPDE toolbox 
%   It tests the second 10 projects of BATCHTEST, with graphical output 
%   Matlab Parallel toolbox is required for a complete test
%   Total runtime should be around 60s, depending on CPU speed.
%   To delete the generated graphs, remove '%' before the  'close all' 
clear;
e(11) = GPE;
e(12) = GPE2;
e(13) = Kubo2;
e(14) = SolitonAv;
e(15) = Nonlinear;
e(16) = GaussianAv2D;
e(17) = GaussianAv4D;
e(18) = GaussianAv6D;
e(19) = GaussianKinetic;
e(20) = SolitonDeriv;
%close all;         %%Deletes all figures if not wanted
fprintf('xSPDE batch test error score = %f \n',prod(e)*10^31);
fprintf('Expected total error score   = 4.294849 \n' );


