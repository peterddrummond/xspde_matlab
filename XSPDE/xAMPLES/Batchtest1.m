function [] = Batchtest1()
%   BATCHTEST1 is a short Matlab function to test the operation of xSPDE 
%   It tests the first 12 projects of the full BATCHTEST script
%   No Matlab parallel toolbox is required
%   Total runtime should be l20-60s, depending on CPU speed.
%   To delete the graphs, remove '%' before the  'close all' command.

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
e(11) = GPE;
e(12) = GPE2;

%close all;         %%Deletes all figures if not wanted

Et = toc(t1);       %%Elapsed time of all simulations

fprintf('\n\nxSPDE error scores, expected results in brackets:\n\n');
e1 = [42.3862,0.000195,0.012153,0.008502,0.000000,67.536143,10.151231,...
    0.007582,2.081024,50.465043,0.000009,0.000935];
for j = 1:12
    fprintf('xSPDE Batchtest # %d error score  = %f (%f)\n',j,e(j),e1(j));
end
fprintf('\nxSPDE Batchtest1, geometric mean of error = %f (0.021498)\n',prod(e)^(1/12));
fprintf('\nxSPDE Batchtest1 elapsed time  = %f seconds (~23) \n',Et);
end