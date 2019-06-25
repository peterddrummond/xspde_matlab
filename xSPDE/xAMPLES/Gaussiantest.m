function [e] = Gaussiantest()
%   e  =  GAUSSIANAV2D() tests xSPDE for a linear Schroedinger equation
%   Tests a (two+one)-dimensional partial differential equation for:
%   (1) Initial gaussian in space
%   (2) Setting da to zero using the zeros function
%   (3) Averaging over the spatial grid
%   (4) Integrating over the spatial grid
%   (5) Integrating over the momentum grid
%   (6) Using transforms, images, transverse, compare, pdimension options
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =        'Gaussian diffraction in 2D';
in.dimension =   3;                                       %%dimension: 1-4
in.initial =     @(w,r) exp(-0.5*(r.x.^2+r.y.^2));        %%initialisation  
in.da =          @(a,~,~) zeros(size(a));                 %%da 
in.observe{1} =  @(a,r) a.*conj(a);                       %%observable 1
in.observe{2} =  @(a,r) a.*conj(a);                       %%observable 1


in.linear =      @(r) 1i*0.05*(r.Dx.^2+r.Dy.^2);          %%laplacian
in.transforms =  {0;[0,1,1]};                               %%transforms
in.olabels =     {'I(x)','I(k)'};%%labels 
in.compare{2} =  @(t,~,~) 1;                              %%comparison 
e  = xspde(in);                                           %%simulation
end                                                       %%end of main 
