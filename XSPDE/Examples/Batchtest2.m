%   BATCHTEST is a Matlab script to test the operation of the xSPDE toolbox 
%   There are currently 20 projects tested, each with graphical output 
%   Matlab Parallel toolbox is required for a complete test
%   Total runtime should be around 60s, depending on CPU speed.
%   To see the generated graphs, add '%' before the  'close all' commands

e(11) = GPE;
e(12) = GPE2;
e(13) = Kubotest;
e(14) = SolitonAv;
e(15) = Nonlinear;
e(16) = GaussianAv2D;
e(17) = GaussianAv4D;
e(18) = GaussianAv6D;
e(19) = GaussianKinetic;
e(20) = SolitonDeriv;
%close all;         %%Deletes all figures if not wanted
fprintf('xSPDE batch test error score = %f \n',prod(e)*10^45);
fprintf('Expected total error score   = 19.354908 \n' );


