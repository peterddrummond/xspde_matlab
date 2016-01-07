function [e] = GaussianAv6D()                                 
%   e  =  GAUSSIANAV6D() tests xSPDE for a linear Schroedinger equation
%   Tests a (six+one)-dimensional partial differential equation for:
%   (1) Initial 6D gaussian in space, using numeric grid axis notation
%   (2) Setting da to zero using the default in 6D
%   (3) Integrating over part of the spatial grid using xint
%   (5) Integrating over part of the the momentum grid
%   (6) Using transforms over part of the spatial grid
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =        'Schroedinger equation in 6 space dimensions';
in.dimension =   7; 
in.points =   [11,11,11,11,11,11,11]; %%dimension: 1-4
in.ranges =   [10,6,6,6,6,6,6];
in.initial =     @(w,r) exp(-0.5*(r.x{1}.^2+r.x{2}.^2+r.x{3}.^2+r.x{4}.^2+r.x{5}.^2+r.x{6}.^2));   
in.observe{1} =  @(a,r) a.*conj(a);                         
in.observe{2} =  @(a,r) xint(a.*conj(a),r.dx.*[0,1,0,1,0,1,1],r);
in.observe{3} =  @(a,r) a.*conj(a);                         
in.observe{4} =  @(a,r) xint(a.*conj(a),[r.dx(1:3),r.dk(4:7)],r);
in.linear =      @(D,r) 1i*0.05*(D{1}.^2+D{2}.^2+D{3}.^2+D{4}.^2+D{5}.^2+D{6}.^2);
in.transforms =  {0,0,[0,0,0,1,1,1,1],[0,0,0,1,1,1,1]};     
in.images =      2;                                         %%number of images
in.imagetype =   2;                                         %%number of images
in.transverse =  1;                                         %%transverse plots
in.raw        =  1;
%in.file       =  'GaussianAv6D.mat';   %%Warning- makes large disk file!
in.olabels =    {'I','\int I dV','I(k)','\int I dK'};
in.compare{1} =  @(t,~) (1+(t/10).^2).^(-3);                %%comparison 
in.compare{2} =  @(t,~) pi^2*(1+(t/10).^2).^(-1);           %%comparison
in.compare{3} =  @(t,~) (1+(t/10).^2).^(-1);                %%comparison 
in.compare{4} =  @(t,~) pi^3+0.*t;                          %%comparison 
in.pdimension = {3,1,1,1,1,1};
e  = xspde(in);                                             %%simulation
end                                                         %%end of main 
