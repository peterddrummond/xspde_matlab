
function [] = Batchtest()
%   BATCHTEST is a Matlab script to test the operation of the xSPDE toolbox 
%   There are currently 25 projects tested, each with graphical output 
%   Matlab Parallel toolbox is required for this test.
%   Runtime: <120s, for R2016a on a Mac, depending on the hardware.
%   To see the generated graphs, add '%' before the  'close all' commands

cd ~               %Use home directory for stored data files

t1=tic;
e{1} = Wiener;           n{1} = 'Wiener';
e{2} = SHO;              n{2} = 'SHO';
e{3} = Kubo;             n{3} = 'Kubo';
e{4} = Soliton;          n{4} = 'Soliton';
e{5} = Gaussian;         n{5} = 'Gaussian';
close all;               %%Deletes all figures if not wanted
e{6} = Planar;           n{6} = 'Planar';
e{7} = Gain;             n{7} = 'Gain';
e{8} = Characteristic;   n{8} = 'Characteristic';
e{9} = Equilibrium;      n{9} = 'Equilibrium';
e{10} = GainTransfer;    n{10} = 'GainTransfer';
close all;               %%Deletes all figures if not wanted
e{11} = GPE;             n{11} = 'GPE';
e{12} = GPE2;            n{12} = 'GPE2';
e{13} = Kubo2;           n{13} = 'Kubo2';
e{14} = SolitonAv;       n{14} = 'SolitonAv';
e{15} = Nonlinear;       n{15} = 'Nonlinear';
close all;               %%Deletes all figures if not wanted
e{16} = GaussianAv2D;    n{16} = 'GaussianAv2D';
e{17} = GaussianAv4D;    n{17} = 'GaussianAv4D';
e{18} = GaussianAv6D;    n{18} = 'GaussianAv6D';
e{19} = GaussianKinetic; n{19} = 'GaussianKinetic';
e{20} = SolitonDerivS;   n{20} = 'SolitonDerivS';
close all;               %%Deletes all figures if not wanted
e{21} = SolitonDerivN;   n{21} = 'SolitonDerivN';
e{22} = SolitonDerivDir; n{22} = 'SolitonDerivDir';
e{23} = GPEvortex2D;     n{23} = 'GPEvortex2D';
e{24} = Quantum;         n{24} = 'Quantum';
e{25} = Kubodefine;      n{25} = 'Kubodefine';
close all;               %%Deletes all figures if not wanted
Et = toc(t1);      
                         %%Print summarized error scores and time

fprintf('\n\nxSPDE3.2 error scores, expected results in brackets (R2019):\n\n');
e1 = [84.428414,0.000115,0.011726,0.008828,0.000000,58.069965,10.674732,...
    0.008131,6.938672,188.033844,0.000015,0.000935,0.417230,0.035262,...
    0.338336,0.000001,0.000055,0.003209,0.000092,0.000247,.000005,...
    0.071056,0.032642,0.211228, 0.315903];
et = 1;
for j = 1:25
    fprintf('Batchtest # %d error  = %f (%f); time = %f (%s)\n',j,e{j}(1),e1(j),e{j}(2),n{j});
    et = et*e{j}(1);
end
%Ex_err = 0.163773;
Ex_err = 0.12;
fprintf('\n Test error score = %f, Expected score = %f - %f\n',et*10^50,Ex_err,2*Ex_err);
if (et*10^50>Ex_err && et*10^50<2*Ex_err) 
 fprintf('\n xSPDE batch test was passed successfully\n');
else
     fprintf('\n xSPDE batch test failed, check installation\n');
end
fprintf('\n Test elapsed time  = %f (<~120) seconds\n',Et);
end