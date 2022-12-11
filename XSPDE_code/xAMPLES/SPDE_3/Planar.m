function [e] = Planar()
%   e  =  PLANAR() tests xSPDE for a linear stochastic PDE
%   Tests a three-dimensional partial stochastic differential equation for:
%   (1) Inputting the ranges in three dimensions
%   (2) Reducing the points in time
%   (3) Defining the noises including noise in Fourier space
%   (4) Inputting more than one stochastic field
%   (5) Using Int and Ave with different numbers of arguments
%   (6) Using the transforms vector to Fourier transform the output fields
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

p.name        =  'Planar noise growth';                       %%name for simulation
p.dimensions  =  3;                                           %%dimension: 1-4  = t,x,y,z
p.fields      =  2;                                           %%field components
p.ranges      =  [1,5,5];                                     %%ranges: t,x,y,z
p.points      =  10;                                          %%steps per plotted point
p.noises      =  [2,2];                                       %%xnoises, knoises per point
p.ensembles   =  [10,2,12];                                   %%samples,ensembles,parallel
p.initial     =  @Initial;                                    %%Initialisation  handle
p.deriv       =  @D_planar;                                   %%Derivative  handle
p.linear      =  @Linear;                                     %%Derivative  handle
p.observe{1}  =  @(a,p) Int(a(1,:).*conj(a(1,:)),p);          %%Observe  handle
p.observe{2}  =  @(a,p) Int(a(2,:).*conj(a(2,:)),p.dk,p);     %%Observe  handle
p.observe{3}  =  @(a,p) real(Ave(a(1,:).*conj(a(2,:)),p));    %%Observe  handle
p.observe{4}  =  @(a,p) Int(a(2,:).*conj(a(2,:)),p);          %%Observe  handle
p.transforms  =  {[0,0,0],[0,1,1],[0,1,1]};
p.olabels{1}  =  '<\int|a_1(x)|^2 d^2x>';                     %%labels
p.olabels{2}  =  '<\int|a_2(k)|^2 d^2k>';                     %%labels
p.olabels{3}  =  '<< a_1(k)a^*_2(k)>>';                       %%labels
p.olabels{4}  =  '<\int|a_2(x)|^2 d^2x>';                     %%labels
p.compare{1}  =  @(p) (1+p.t)*p.nspace;
p.compare{2}  =  @(p) (1+p.t)*p.nspace;
p.compare{3}  =  @(p) 0.0;
p.pdimension  =  {1,1,1,1};                                   %%maximum plot dimension
e             =  xspde(p);                                    %%Stochasic program
end                                                           %%end of main function

							      %%XSPDE user functions

function a0 = Initial(v,~)                                    %%Initialises fields
a0(1,:)  =  (v(1,:)+1i*v(2,:))/sqrt(2);                       %%lattice vectors
a0(2,:)  =  (v(3,:)+1i*v(4,:))/sqrt(2);                       %%lattice vectors
end                                                           %%end initialise fields

function da = D_planar(~,w,~)                                 %%Derivatives
da(1,:)  =  (w(1,:)+1i*w(2,:))/sqrt(2);                       %%complex x-noise equation
da(2,:)  =  (w(3,:)+1i*w(4,:))/sqrt(2);                       %%complex k-noise equation
end                                                           %%end local derivatives

function L = Linear(p)                                        %%Linear coefficient
lap     =  p.Dx.^2+p.Dy.^2;                                   %%Laplacian
L(1,:)  =  1i*0.5*lap(:);                                     %%Linear coefficient(1)
L(2,:)  =  1i*0.5*lap(:);                                     %%Linear coefficient(2)
end                                                           %end linear coefficient
