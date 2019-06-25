function e = SHO()                            
%   e  =  SHO() tests xSPDE for a simple harmonic oscillator
%   Tests a one-variable differential equation for: 
%   (1) Vector calculations
%   (2) Exactly soluble DE case
%   (3) Default range, points, dimension
%   (4) Adding an initial value
%   (5) Default observe functions, olabels
%   (6) Comparisons, scatters,parametric plots, sequential plots
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

 in.name =      'Simple harmonic oscillator';    %%Name 
 in.initial =   @(w,~) [1,1,1,;0,0,0];                    %%Initialisation
 in.fields =    2;
 in.steps =    20;
 in.ensembles = 3;
 in.olabels   = {'x'};
 in.da =        @(a,w,r) [-a(2,:);a(1,:)];      %%Derivative function   
 in.compare =   {@(t,in) cos(t)};               %%Compare  handle
 in.observe{1} = @(a,r) a(1,:);
 in2 = in;
 in2.da =        @(a,w,r) [-a(2,:);a(1,:)]+0.005*w;%%Derivative function
 in2.compare{1} = [];
 in2.observe{2} = @(a,r) a(2,:);
 in2.observe{3} = @(a,r) a(1,:);
 in2.scatters =  {1,3,3};
 in2.olabels   = {'x';'y';'x'};
 in2.parametric{3} = [2,1];
 e =xspde({in,in2});
end                                             %%end of main function
