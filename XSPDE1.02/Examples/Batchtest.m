%%BATCHTEST script for XSPDE1.0
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
e(11) = GPE;
e(12) = GPE2;
e(13) = Kubotest;
e(14) = SolitonAv;
e(15) = Nonlinear;
e(16) = GaussianAv1D;
e(17) = GaussianAv2D;
e(18) = GaussianAv3D;
close all;         %%Deletes all figures if not wanted
fprintf('Batch test error score = %f \n',prod(e)*10^40);
%%Expected answer:  Batch test error score = 6.482868 

