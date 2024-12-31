function [] = Batchtest()
%   BATCHTEST is a  script to test the operation of the xSPDE toolbox 
%   There are currently 42 projects tested, each with graphical output 
%   Matlab Parallel toolbox is required for this test.
%   If not available, replace parallel ensembles by serial ensembles
%   Runtime: 20-100s for Matlab R2024, depending on cores available
%            <1200s for Octave v6.2
%   To see the generated graphs, add '%' before the  'close all' commands
%   If needed, you can also include: 
%cd ~               %Use home directory for stored data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=tic;
e{1} = Wiener;           n{1} = 'Wiener';
e{2} = SHO;              n{2} = 'SHO';
e{3} = Wienercell;       n{3} = 'Wienercell';
e{4} = Soliton;          n{4} = 'Soliton';
e{5} = Gaussian3D;       n{5} = 'Gaussian3D';
close all;               %%Deletes all figures if not wanted
e{6} = Planarstep;       n{6} = 'Planarstep';
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
e{26} = Wienertest;      n{26} = 'Wienertest';
e{27} = Wienertest2;     n{27} = 'Wienertest2';
e{28} = Weight;          n{28} = 'Weight';
e{29} = Subharmonic;     n{29} = 'Subharmonic';
e{30} = Wienerprob;      n{30} = 'Wienerprob' ;
close all;               %%Deletes all figures if not wanted
e{31} = Catenoid;        n{31} = 'Catenoid' ;
e{32} = Algorithms;      n{32} = 'Algorithms';
e{33} = Boundaries2;     n{33} = 'Boundaries2';
e{34} = Boundaries3;     n{34} = 'Boundaries3';
e{35} = Sphere10;        n{35} = 'Sphere10';
e{36} = SSElin;          n{36} = 'SSElin';
e{37} = SSEnonlin2spr;   n{37} = 'SSEnonlin2spr';
e{38} = Masterlin2;      n{38} = 'Masterlin2';
e{39} = SSElin2spr;      n{39} = 'SSElin2spr';
e{40} = SSElin2j;        n{40} = 'SSElin2j';
e{41} = SSElin2jsp;      n{41} = 'SSElin2jsp';
e{42} = Peregrine;       n{42} = 'Peregrine';

close all;               %%Deletes all figures if not wanted
Et = toc(t1);      
                         %%Print summarized error scores and time

fprintf('\n\nxSPDE4.2 Batchtest RMS errors, expected in brackets:\n\n');
fprintf(['Result: (Expected)     Total        Step      ',...
'Sampling  Diff.     Chi-squ/k Timing    Name\n\n']);
e1 = [9.720820e-03,1.640870e-05,4.119120e-03,5.047219e-03,1.682525e-11,...
      1.360118e-02,2.031477e-03,9.654854e-05,2.436846e-02,9.409008e-02,...
      4.069716e-06,6.585455e-03,4.266752e-02,1.410144e-02,1.455296e-02,...
      4.813746e-07,2.444070e-07,6.243503e-05,4.049843e-06,8.570537e-08,...
      2.134464e-06,6.292213e-04,1.772071e-03,1.454942e-02,1.171775e-02,...
      1.223001e-02,1.522842e-03,1.124233e-02,1.244600e-03,3.763128e-03,...
      8.133844e-03,5.328096e-04,2.185585e-05,1.105256e-03,3.069677e-03,...
      3.296567e-03,1.911476e-03,4.236356e-08,3.030441e-03,1.798514e-02,...
      1.327488e-02,4.588574e-04];
et = zeros(1,6);
nt = zeros(1,6);
Exp = 1;
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
  nt = (e{j} > 1.e-15) + nt;
  et(5:6) = et(5:6)+e{j}(5:6);
  Exp = Exp*e1(j);
end
et(1) = et(1)^(1/nt(1));
et(2:5) = et(2:5)./nt(2:5);
et(2:4) = sqrt(et(2:4));
Exp = Exp^(1/nt(1));
fprintf('\nBatchtest errors vs expected:\n\n');
fprintf('Total (GM)     = %d, tests = %d, Expect = %d\n',et(1),nt(1),Exp);
fprintf('Step (RMS)     = %d, tests = %d\n',et(2),nt(2));
fprintf('Sampling (RMS) = %d, tests = %d\n',et(3),nt(3));
fprintf('Diff. (RMS)    = %d, tests = %d\n',et(4),nt(4));
fprintf('Chi-square/k   = %d, tests = %d\n',et(5),nt(5));
fprintf('\nTotal time     = %.3g seconds\n',et(6));
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