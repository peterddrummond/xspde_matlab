function [e] = Algorithms()
%   e  =  Alorithms() tests xSPDE for a nonlinear Schroedinger equation
%   using 7 different builtin xSPDE algorithms
%   Tests a one-dimensional stochastic partial differential equation using:
%   (1) MP
%   (2) Implicit
%   (3) MPadapt
%   (4) RK2
%   (5) RK4
%   (6) Euler
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License

p.name       =  'Algorithms:MP';                        %%Graph title
p.dimensions =  2;                                      %%Space-time dimn
p.ensembles  =  [6,6];                                 %%Ensembles
p.points     =  [31,31];                                %%Points
p.ranges     =  [5,15];                                 %%Ranges
p.steps      =  5;
p.initial    =  @(w,p)   2*sech(p.x);                   %%Initialisation
p.transfer   =  @(~,w,p) 2*sech(p.x);                   %%Initialisation
p.deriv      =  @(a,w,~) 1i*(.005*w+a.*conj(a)/4).*a;   %%Nonlinear deriv
p.linear     =  @(p)     0.5*1i*(p.Dx.^2-1.0);          %%Linear terms
p.method     =  @MP;
p1           =  p;
p2           =  p;
p2.steps     =  40;
p2.name      =  'Algorithms:Implicit';               
p2.method    =  @Implicit;
p3           =  p;
p3.method    =  @MPadapt;
p3.name      =  'Algorithms:MPadapt';                  
p4           =  p;
p4.method    =  @RK2;
p4.name      =  'Algorithms:RK2';                     
p5           =  p;
p5.method    =  @RK4;
p5.name      =  'Algorithms:RK4';                    
p6           =  p;
p6.method    =  @Euler;
p6.steps     =  40;
p6.name      =  'Algorithms:Euler';  
e            =  xspde(p1,p2,p3,p4,p5,p6);             %%main program
end                                                     %%end of function
