function [e] = GaussianAv4D()
%   e  =  GAUSSIANAV4D() tests xSPDE for a linear Schroedinger equation
%   Tests a (four+one)-dimensional partial differential equation for:
%   (1) Initial 4D gaussian in space, using numeric grid axis notation
%   (2) Setting da to zero using the default
%   (3) Averaging over the spatial grid using one argument for xave
%   (4) Integrating over the spatial grid using two arguments for xint
%   (5) Averaging and integrating over the momentum grid
%   (6) Using transforms, images, transverse, compare, pdimension in 4+1 D
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

cd ~
in.name =        '4+1D Schroedinger equation';
in.dimension =   5;
in.noises =      2;
in.points =      [7,15,15,15,15];                         
in.initial =     @(w,r) exp(-0.5*(r.x{2}.^2+r.x{3}.^2+r.x{4}.^2+r.x{5}.^2));  
in.observe{1} =  @(a,r) a.*conj(a);                         %%observable 1
in.observe{2} =  @(a,r) xave(a.*conj(a));                   %%observable 2
in.observe{3} =  @(a,r) xint(a.*conj(a),r);                 %%observable 3
in.observe{4} =  @(a,r) a.*conj(a);                         %%observable 4
in.observe{5} =  @(a,r) xave(a.*conj(a));                   %%observable 5
in.observe{6} =  @(a,r) xint(a.*conj(a),r.dk,r);            %%observable 6
in.linear =      @(r) 1i*0.05*(r.D{2}.^2+r.D{3}.^2+r.D{4}.^2+r.D{5}.^2);
in.transforms =  {0,0,0,[0,1,1,1,1],[0,1,1,1,1],[0,1,1,1,1]};     
in.images =      2;                                         %%number of images
in.transverse =  {2,0,0,2};                                 %%transverse plots
in.olabels =    {'I','<I>','\int I dV','I(k)','<I(k)>','\int I dK'};%%labels 
in.compare{1} =  @(t,~) (1+(t/10).^2).^(-2);                %%comparison 
in.compare{4} =  @(t,~) 1+0.*t;                             %%comparison 
in.pdimension = {3,1,1,2,1,1};
%in.file       = 'GaussianAv4D.h5';%%Warning - makes large disk file
e =             xspde(in); 
end                                                         %%end of main 
