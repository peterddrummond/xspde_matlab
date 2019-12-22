function [e] = GPE2()                                    
%   e  =  GPE2() tests xSPDE for  a GPE in a harmonic well
%   Tests a two-dimensional PDE with two fields
%   (1) Modifying points and ranges in 2 dimensions
%   (2) Increasing steps between points
%   (3) Setting the initial fields for a two component field
%   (4) Using an auxiliary function of space to define a potential
%   (5) Defining a linear response for two fields
%   (6) Defining function handles separately from the derivative
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =       'GPE2';
in.dimension =  2; 
in.fields =     2;  
in.points =     [101,64];
in.ranges =     [20,20];
in.steps  =     4;
in.initial =    @(w,r)          [sech(r.x+5);sech(r.x-5)];
V          =    @(r)            0.1*r.x.^2;
rho        =    @(a)            conj(a(1,:)).*a(1,:)+conj(a(2,:)).*a(2,:);
da1        =    @(a,r)           -1i*a(1,:).*(V(r)+rho(a));
da2        =    @(a,r)           -1i*a(2,:).*(V(r)+rho(a));
in.da =         @(a,xi,r)       [da1(a,r);da2(a,r)];
in.linear =     @(r)          [0.5*1i*(r.Dx.^2);0.5*1i*(r.Dx.^2)]; 
in.observe{1} = @(a,r) a(1,:).*conj(a(1,:));
in.observe{2} = @(a,r) a(2,:).*conj(a(2,:));
in.images     = {2};                                       
in.imagetype  = {2};  
in.olabels =    {'|a_1|^2','|a_2|^2'};
e  =            xspde(in);                             %%main program
end                                                    %%end of function
