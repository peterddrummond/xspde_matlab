function [e] = GPE2()
%   e  =  GPE2() tests xSPDE for  a GPE in a harmonic well
%   Tests a two-dimensional PDE with two fields
%   (1) Modifying points and ranges in 2 dimensions
%   (2) Increasing steps between points
%   (3) Initial fields for a two components, and vector parallel ensembles
%   (4) Using an auxiliary function of space to define a potential
%   (5) Defining a linear response for two fields
%   (6) Defining function handles separately from the derivative
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name         =  'GPE2';
p.dimensions   =  2;
p.fields       =  2;
p.order        =  2;
p.points       =  [101,64];
p.ranges       =  [20,20];
p.ensembles    =  6;
p.steps        =  4;
p.initial      =  @(w,p)     [sech(p.x+5);sech(p.x-5)];
V              =  @(p)       0.1*p.x.^2;
rho            =  @(a,p)     xshape(conj(a(1,:)).*a(1,:)...
                             +conj(a(2,:)).*a(2,:),1,0,p);
p.deriv        =  @(a,w,p)  -1i*a.*([V(p);2*V(p)]+rho(a,p))+0.01*w;
p.linear       =  @(p)       [0.5*1i*(p.Dx.^2);0.5*1i*(p.Dx.^2)];
p.observe{1}   =  @(a,p)     a(1,:).*conj(a(1,:));
p.observe{2}   =  @(a,p)     a(2,:).*conj(a(2,:));
p.images       =  {2};
p.imagetype    =  {2};
p.olabels      =  {'|a_1|^2','|a_2|^2'};
e              =  xspde(p);                           %%main program
end                                                   %%end of function
