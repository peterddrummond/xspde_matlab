function [e] = GPEvortex2D()
%   e  =  GPEvortex2D() tests xSPDE for a lattice of superfluid vortices
%   Tests a 2D GPE with vortices plus rotation
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name         =  'GPEvortex2D';
p.dimensions   =  3;
p.fields       =  1;
p.points       =  [50,40,40];
p.ranges       =  [15,16,16];
p.steps        =  15;
g              =  200;
om             =  0.6;
L              =  @(a,p)       1i*(p.x.*D1(a,3,p)-p.y.*D1(a,2,p));
V              =  @(p)         0.35*(p.x.^2+p.y.^2);
p.initial      =  @(v,p)       0.1*exp(-V(p));
rho            =  @(a)         g*conj(a).*a;
p.deriv        =  @normda;
p.da1          =  @(a,w,p)     -a.*(V(p)+rho(a))+om*L(a,p);
p.linear       =  @(p)         0.5*(p.Dx.^2+p.Dy.^2);
p.observe{1}   =  @(a,p) a(1,:).*conj(a(1,:));
p.observe{2}   =  @(a,p) a(1,:).*conj(a(1,:));
p.images       =  {2,2};
p.imagetype    =  {1,2};
p.olabels      =  {'|a|^2','|a|^2'};
e              =  xspde(p);                             %%main program

function b = normda(a,w,p)
%   b  =  NORMDA(a,z,r) is a normalized derivative
%   Takes a derivative and returns a normalized step
b     =  a+p.da1(a,w,p)*p.dtr;
norm  =  sqrt(Int(abs(b).^2,p.dx,p));
b     = (b./norm-a)/p.dtr;
end
end                                                    %%end of function
