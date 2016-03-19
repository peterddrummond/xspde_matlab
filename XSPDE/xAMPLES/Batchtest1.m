%   BATCHTEST1 is a short Matlab script to test the operation of xSPDE 
%   It tests the first 10 projects of the full BATCHTEST script
%   Matlab Parallel toolbox is required for a complete test
%   Total runtime should be less than 60s, depending on CPU speed.
%   To delete the graphs, remove '%' before the  'close all' commands

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
%close all;         %%Deletes all figures if not wanted
fprintf('xSPDE Batch1 test error score = %f \n',prod(e)*10^11);
fprintf('Expected total error score   = 0.768651 \n' );


