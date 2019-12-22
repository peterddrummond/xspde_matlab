function [e] = GPEvortex2D1()                                    
%   e  =  GPEvortex2D() tests xSPDE for a lattice of superfluid vortices
%   Tests a 2D GPE with vortices plus rotation 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =       'GPEvortex2D1';
in.dimension =  3; 
in.fields =     1;  
in.points =     [50,63,63];
in.ranges =     [15,20,20];
in.steps =      15;
g          =    200;
om =            0.6;
%L          =    @(a,r)       1i*(r.x.*xd(a,r.Dy,r)-r.y.*xd(a,r.Dx,r));
V          =    @(r)         0.35*(r.x.^2+r.y.^2);
in.initial =    @(w,r)       0.1*exp(-V(r));
rho        =    @(a)         g*conj(a(1,:)).*a(1,:);
in.da     =     @normda;
in.da1     =    @(a,xi,r)     -a(1,:).*(V(r)+rho(a))+om*L(a,r);
in.linear =     @(r)        0.5*(r.Dx.^2+r.Dy.^2); 
in.observe{1} = @(a,r) a(1,:).*conj(a(1,:));
in.observe{2} = @(a,r) a(1,:).*conj(a(1,:));
in.images    =  {2,2}; 
in.imagetype  = {1,2}; 
in.olabels =    {'|a_1|^2','|a_1|^2'};
e  =            xspde(in);                             %%main program

function b = normda(a,~,r)
%   b  =  NORMDA(a,~,r) is a normalized derivative
%   Takes a derivative and returns a normalized step

da1 =-a(1,:).*(0.35*(r.x.^2+r.y.^2)+rho(a));
da1 =da1+om*1i*(r.x.*xd(a,r.Dy,r)-r.y.*xd(a,r.Dx,r));
b = a+da1*r.dtr;
norm = sqrt(xint(abs(b).^2,r.dx, r));
b = (b./norm-a)/r.dtr;
end
end                                                    %%end of function



