%   BATCHTEST is a Matlab script to test the operation of the xSPDE toolbox 
%   There are currently 21 projects tested, each with graphical output 
%   Matlab Parallel toolbox is required for a complete test
%   Total runtime should be around 60s, depending on CPU speed.
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
e(13) = Kubo2;
e(14) = SolitonAv;
e(15) = Nonlinear;
e(16) = GaussianAv2D;
e(17) = GaussianAv4D;
e(18) = GaussianAv6D;
e(19) = GaussianKinetic;
e(20) = SolitonDeriv;
e(21) = GPEvortex2D;
close all;         %%Deletes all figures if not wanted

fprintf('xSPDE error scores, expected reults in brackets:\n');
e1 = [0.185807,0.034223,0.012153,0.037387,0.000000,74.251926,19.535136,...
    0.008131,4.858198,50.560601,0.000009,0.076435,0.237978,0.065598,...
    0.338336,0.000001,0.000001,0.017825,0.000092,0.000123,0.001375];
for j = 1:21
    fprintf('xSPDE test # %d error score  = %f (%f)\n',j,e(j),e1(j));
end
fprintf('xSPDE batch test error score = %f (11.549549)\n',prod(e)*10^45);