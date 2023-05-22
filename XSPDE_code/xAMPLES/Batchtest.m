function [] = Batchtest()
%   BATCHTEST is a  script to test the operation of the xSPDE toolbox 
%   There are currently 35 projects tested, each with graphical output 
%   Matlab Parallel toolbox is required for this test.
%   If not available, replace parallel ensembles by serial ensembles
%   Runtime: <120s for Matlab R2020, depending on cores available
%            <1200s for Octave v6.2
%   To see the generated graphs, add '%' before the  'close all' commands
%   If needed, you can also include: 
%cd ~               %Use home directory for stored data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=tic;
e{1} = Wiener;           n{1} = 'Wiener';
e{2} = SHO;              n{2} = 'SHO';
e{3} = Kubo;             n{3} = 'Kubo';
e{4} = Soliton;          n{4} = 'Soliton';
e{5} = Gaussian3D;       n{5} = 'Gaussian3D';
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
e{26} = Wienertest;      n{26} = 'Wienertest';
e{27} = Wienertest2;     n{27} = 'Wienertest2';
e{28} = Weight;          n{28} = 'Weight';
e{29} = Subharmonic;     n{29} = 'Subharmonic';
e{30} = Wienerprob;      n{30} = 'Wienerprob' ;
e{31} = Catenoid;        n{31} = 'Catenoid' ;
e{32} = Algorithms;      n{32} = 'Algorithms';
e{33} = Boundaries2;     n{33} = 'Boundaries2';
e{34} = Boundaries3;     n{34} = 'Boundaries3';
e{35} = Sphere10;        n{35} = 'Sphere10';

close all;               %%Deletes all figures if not wanted
Et = toc(t1);      
                         %%Print summarized error scores and time

fprintf('\n\nxSPDE3.44 Batchtest RMS errors, expected in brackets:\n\n');
fprintf(['Result: (Expected)     Total        Step      ',...
'Sampling  Diff.     Chi-squ/k Timing    Name\n\n']);
e1 = [1.128155e-02,1.342241e-05,4.724121e-03,3.522462e-03,1.682525e-11,...
    1.544402e-03,2.008696e-03,9.287161e-05,1.603171e-02,9.404619e-02,...
    2.998764e-06,2.017065e-02,4.866021e-02,1.439964e-02,1.455296e-02,...
    4.813751e-07,2.444071e-07,6.243846e-05,4.049878e-06, 8.522040e-08,...
    2.134473e-06,8.260400e-05,1.773850e-03,1.381911e-02,9.686027e-03,...
    1.209275e-02,5.236364e-03,5.997868e-03,1.975394e-03,4.765677e-03,...
    7.915032e-03,6.071103e-04,2.186477e-05,1.105256e-03, 2.779772e-03];
et = zeros(1,6);
nt = zeros(1,6);
et(1) = 1;
for j = 1:length(e)
  if e{j}(1)>e1(j)/1.001 && e{j}(1)<e1(j)*1.001 
    fprintf('# %2d OK (%.6e) %.6e %.3e %.3e %.3e %.3e %.3e (%s)\n',...
        j,e1(j),e{j}(1:6),n{j});
  else
    fprintf('# %2d ?? (%.6e) %.6e %.3e %.3e %.3e %.3e %.3e (%s)\n',...
        j,e1(j),e{j}(1:6),n{j});
  end
  et(1) = et(1)*e{j}(1);
  et(2:4) = et(2:4)+e{j}(2:4).^2;
  nt = (e{j} > 1.e-99) + nt;
  et(5:6) = et(5:6)+e{j}(5:6);
end
et(1) = et(1)^(1/nt(1));
et(2:5) = et(2:5)./nt(2:5);
et(1:4) = sqrt(et(1:4));
Exp = 1.939832e-02;
fprintf('\nBatchtest errors vs expected:\n');
fprintf('\nTotal (GM)       = %d, tests = %d, Expected = %d\n',et(1),nt(1),Exp);
fprintf('Step (RMS)       = %d, tests = %d\n',et(2),nt(2));
fprintf('Sampling (RMS)   = %d, tests = %d\n',et(3),nt(3));
fprintf('Difference (RMS) = %d, tests = %d\n',et(4),nt(4));
fprintf('Av. chi-square/k = %d, tests = %d\n',et(5),nt(5));
fprintf('\nSimulation time  = %.3g seconds\n',et(6));
if (et(1)>Exp/1.001 && et(1)<1.001*Exp) 
 fprintf('\nTotal error within 0.1%% of expected error\n');
else
  if (et(1)>Exp/1.05 && et(1)<1.05*Exp) 
    fprintf('\nTotal error within 5%% of expected value\n');
    fprintf('If Octave, errors may differ from Matlab error estimates \n');
  else
    fprintf('\nBatchtest failed, check installation\n');
  end
end
fprintf('\nTotal elapsed time  = %.3g (<~120) seconds\n',Et);
end