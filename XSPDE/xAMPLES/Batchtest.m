
function [] = Batchtest()
%   BATCHTEST is a Matlab script to test the operation of the xSPDE toolbox 
%   There are currently 22 projects tested, each with graphical output 
%   Matlab Parallel toolbox is required for this test.
%   Runtime: <400s, for R2016a on a Mac, depending on the hardware.
%   To see the generated graphs, add '%' before the  'close all' commands

cd ~               %Use home directory for stored data files

t1=tic;
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
e(22) = Quantum;
close all;         %%Deletes all figures if not wanted
Et = toc(t1);      

                   %%Print summarized error scores and time

fprintf('\n\nxSPDE2.4 error scores, expected reults in brackets:\n\n');
e1 = [84.428414,0.000195,0.012153,0.008502,0.000000,58.626392,10.151231,...
    0.008131,6.705190,50.560601,0.000007,0.000935,0.237978,0.035262,...
    0.338336,0.000001,0.000001,0.016355,0.000092,0.002250,0.024137,0.052807];
for j = 1:22
    fprintf('Batchtest # %d error score  = %f (%f)\n',j,e(j),e1(j));
end
fprintf('\n Test error score = %f (0.290921)\n',prod(e)*10^45);
fprintf('\n Test elapsed time  = %f (<~400) seconds\n',Et);
end