function [e] = Planarstep()
%   e  =  PLANAR() tests xSPDE for a linear stochastic PDE
%   Tests a three-dimensional partial stochastic differential equation for:
%   (1) Inputting the ranges in three dimensions
%   (2) Reducing the points in time to 10, and using a steps input in space
%   (3) Defining the noises including noise in Fourier space
%   (4) Inputting more than one stochastic field
%   (5) Using Int and Ave with different numbers of arguments
%   (6) Using the transforms vector to Fourier transform the output fields
%   Licensed by Peter D. Drummond, (2024) - see License

p.name        =  'Planar noise growth';          %%name for simulation
p.dimensions  =  3;                              %%dimension: 1-4 = t,x,y,z
p.fields      =  2;                              %%field components
p.ranges      =  [1,5,5];                        %%ranges: t,x,y
p.steps       =  [1,1,2];                        %%steps per plotted point
p.points      =  10;                             %%time points
p.noises      =  2;                              %%xnoises, per point
p.knoises     =  2;                              %% knoises per point
p.inrandoms   =  2;                              %%xnoises, per point
p.krandoms    =  2;                              %% knoises per point
p.ensembles   =  [10,2,12];                      %%samples,serial,parallel                   
p.initial     =  @Initial;                       %%Initialisation  handle
p.deriv       =  @D_planar;                      %%Derivative  handle
p.linear      =  @(p) 0.5*1i*(p.Dx.^2+p.Dy.^2);  %%Laplacian
p.nfilter     =  @(w,p) w.*exp(-p.kx.^2/10);
p.observe{1}  =  @(a,p) Int(a(1,:).*conj(a(1,:)),p);       %%Observe handle
p.observe{2}  =  @(a,p) Int(a(1,:).*conj(a(1,:)),p.dk,p);  %%Observe handle
p.observe{3}  =  @(a,p) real(Ave(a(1,:).*conj(a(2,:)),p)); %%Observe handle
p.observe{4}  =  @(a,p) Int(a(2,:).*conj(a(2,:)),p);       %%Observe handle
p.observe{5}  =  @(a,p) a(2,:).*conj(a(2,:));              %%Observe handle
p.transforms  =  {[0,0,0],[0,1,1],[0,1,1],[0,0,0],[0,1,1]};
p.olabels{1}  =  '<\int|a_1 (x)|^2 d^2 x>';      %%labels
p.olabels{2}  =  '<\int|a_2 (k)|^2 d^2 k>';      %%labels
p.olabels{3}  =  '<< a_1 (k)a^*_2 (k)>>';        %%labels
p.olabels{4}  =  '<\int|a_2 (x)|^2 d^2 x>';      %%labels
p.olabels{5}  =  '<|a_1 (k)|^2>';                %%labels
p.compare{1}  =  @(p) (1+p.t)*p.nspace;
p.compare{2}  =  @(p) (1+p.t)*p.nspace;
p.compare{5}  =  @(p) (1+p.t.*exp(-p.kx.^2/5))*p.v/(2*pi)^2;
p.compare{3}  =  @(p) 0.0;
p.pdimension  =  {1,1,1,1};                      %%maximum plot dimension
p.transverse{5}= 1;
p.images{5}  = 1;
e             =  xspde(p);                       %%Stochasic program
end                                              %%end of main function

							      %%XSPDE user functions

function a0 = Initial(u,v,~)                     %%Initialises fields
a0(1,:,:)  =  (u(1,:,:)+1i*u(2,:,:))/sqrt(2);    %%lattice vectors
a0(2,:,:)  =  (v(1,:,:)+1i*v(2,:,:))/sqrt(2);    %%lattice vectors
end                                              %%end initialise fields

function da = D_planar(~,w,z,~)                  %%Derivatives
da(1,:)  =  (w(1,:)+1i*w(2,:))/sqrt(2);          %%complex x-noise equation
da(2,:)  =  (z(1,:)+1i*z(2,:))/sqrt(2);          %%complex k-noise equation
end                                              %%end local derivatives