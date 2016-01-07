function e = SHO()                            
%   e  =  SHO() tests xSPDE for a simple harmonic oscillator
%   Tests a one-variable differential equation for: 
%   (1) Complex variable calculations
%   (2) Exactly soluble DE case
%   (3) Default range, points,fields, dimension
%   (4) Adding an initial value
%   (5) Extended print
%   (6) Default observe functions, olabels and comparisons
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

 in.name =      'Simple harmonic oscillator';   %%Name 
 in.initial =   @(~,~) 1;                       %%Initialisation
 in.print =     2;
 in.da =        @(a,~,r) 1i*a;                  %%Derivative function   
 in.compare =   {@(t,in) cos(t)};               %%Compare  handle
 e =xspde(in);
end                                             %%end of main function
