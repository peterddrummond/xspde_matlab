function [e] = Subharmonic()
%   [e]  =  Subharmonic() tests xSPDE for a subharmonic quantum oscillator
%   uses a positive-P quasi-distribution mapping
%   plots 6 scatter plot trajectories, demonstrating symmetry breaking
%   Shows quantum tunneling in open systems
%   Mean particle number is about 0.4
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name         =  'Subharmonic';
p.ensembles    =  [10,1];                       %%samples,ensembles,parallel
p.initial      =  @(w,~)    0+0*w ;             %%Initialisation  handle
p.steps        =  200;
p.ranges       =  25;
p.fields       =  2;
p.C            =  1;
p.L            =  4;
p.deriv        =  @Subharmonicda;               %%Derivative  handle
p.observe{1}   =  @(a,~) (a(1,:)+a(2,:))/2;
p.observe{2}   =  @(a,~) a(1,:).*a(2,:);
p.scatters     =  {6};
p.olabels      =  {'x quadrature','particle number'};
p.xlabels      =  {'\tau'};
e              =  xspde(p);                     %%Stochasic program
end

function da = Subharmonicda(a,w,p)
a2   =  (1-a.*a);
a2s  =  sqrt(a2);
b    =  flip(a,1);
da   =  -(p.C+0.5)*a+p.L*b.*a2+a2s.*w;
end
