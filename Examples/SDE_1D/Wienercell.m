function [e]       =  Wienercell()
%   [e]  =  WIENER() tests xSPDE for a basic Wiener process
%   Tests a matrix additive SDE equation for:
%   (1) Initial uniform noises
%   (2) Using two cells of vector noises
%   (3) Reshaping the output of observe
%   (4) Inputting an inline derivative with two noise cells
%   (5) Computing the variance of one element
%   (6) Computing the correlation of two element
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

 p.name       = 'Wienercell';                %%name of simulation
 p.method     = @Euler;                      %%first method
 p.urandoms   = [2,2];
 p.fields     = [2,2];
 p.noises     = {[2,1],[1,2]};
 p.deriv      = @(a,u,v,p) u+v;
 p.initial    = @(u,p) 4*u;
 p.ensembles  = [1000,4,10];                 %%ensembles used
 p.observe{1} = @(a,p) a.^2; 
 p.observe{2} = @(a,p) a(1,1,:).*a(1,2,:); 
 p.compare{1} = @(p) 16/3 + 2*p.t;           %%Comparison
 p.compare{2} = @(p) 4 + p.t;                %%Comparison
 p.olabels{1} = '< a_{11}^2 >';              %%labels
 p.olabels{2} = '< a_{11}a_{12} >';          %%labels
 e            = xspde(p);                    %%Runs xspde simulation
end                                          %%end of main function