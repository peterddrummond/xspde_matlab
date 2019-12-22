function [e] = Gaussian()  
%   e  =  GAUSSIAN() tests xSPDE for a linear Schroedinger equation
%   Tests a four-dimensional partial differential equation for:
%   (1) Inputting the x,y,z-grid into the initial condition
%   (2) Default derivative function set to zero in a PDE
%   (3) Setting up the 3D linear response fincluding a linear loss term
%   (4) Inputting compare function
%   (5) Graphing results with transverse plots
%   (6) Graphing results with movie images
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =       'Gaussian diffraction with loss';
in.dimension =  4;                                         %%dimension: 1-4
in.initial =    @(w,r) exp(-0.5*(r.x.^2+r.y.^2+r.z.^2));   %%initialisation  
in.observe{1} = @(a,~) a.*conj(a);                         %%x-observable
in.linear =     @(r) -0.5+1i*0.05*(r.Dx.^2+r.Dy.^2+r.Dz.^2);%%laplacian
in.images =     {2};                                       %%number of images
in.transverse = {2};                                       %%transverse plots
in.olabels =    {'|a(r)|^2'};                              %%labels
in.compare =    @(t,~) (1+(t/10).^2).^(-3/2).*exp(-t);     %%comparison fn 
e  = xspde(in);                                            %%simulation
end                                                        %%end of main 
