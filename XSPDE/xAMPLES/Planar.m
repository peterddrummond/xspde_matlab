function [e] = Planar()                        
%   e  =  PLANAR() tests xSPDE for a linear stochastic PDE
%   Tests a three-dimensional partial stochastic differential equation for:
%   (1) Inputting the ranges in three dimensions
%   (2) Inputting the steps for greater accuracy
%   (3) Defining the noises including noise in Fourier space
%   (4) Inputting more than one stochastic field
%   (5) Testing xint and xave with different numbers of arguments
%   (6) Using the transforms vector to Fourier transform the output fields
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =       'Planar noise growth';         %%name for simulation
in.dimension =  3;                             %%dimension: 1-4 = t,x,y,z
in.fields =     2;                             %%field components
in.ranges =     [1,5,5];                       %%ranges: t,x,y,z  
in.steps =      2;                             %%steps per plotted point
in.step  =      @xMP;
in.noises =     [4,2];                         %%xnoises, knoises per point
in.ensembles =  [10,10,1];                     %%samples,ensembles,parallel
in.initial =    @Initial;                      %%Initialisation  handle
in.da  =        @D_planar;                     %%Derivative  handle
in.linear =     @Linear;                       %%Derivative  handle
in.observe{1} = @(a,r) xint(a(1,:).*conj(a(1,:)),r);%%Observe  handle
in.observe{2} = @(a,r) xint(a(2,:).*conj(a(2,:)),r.dk,r); %%Observe  handle
in.observe{3} = @(a,r) xave(a(1,:).*conj(a(2,:)),r);%%Observe  handle
in.transforms = {[0,0,0],[0,1,1],[0,1,1]};
in.olabels{1} = '<\int|a_1(x)|^2 d^2x>';       %%labels  
in.olabels{2} = '<\int|a_2(k)|^2 d^2k>';       %%labels
in.olabels{3} = '<< a_1(k)a^*_2(k)>>';         %%labels 
in.compare{1} = @(t,in) [1+t]*in.nspace;
in.compare{2} = @(t,in) [1+t]*in.nspace;
in.compare{3} = @(t,in) 0*t;
in.pdimension = {1,1,1,1};                     %%maximum plot dimension
e  =  xspde(in);                               %%Stochasic program
end                                            %%end of main function

%%XSPDE user functions

function a0 = Initial(w,r)                     %%Initialises fields
a0(1,:) = (w(5,:)+1i*w(6,:))/sqrt(2);          %%lattice vectors
a0(2,:) = (w(3,:)+1i*w(4,:))/sqrt(2);          %%lattice vectors
end                                            %%end initialise fields

function da  =  D_planar(a,z,r)                %%Derivatives
da(1,:) = (z(5,:)+1i*z(6,:))/sqrt(2);          %%complex k-noise equation
da(2,:) = (z(3,:)+1i*z(4,:))/sqrt(2);          %%complex x-noise equation
end                                            %%end local derivatives

function L = Linear(r)                         %%Linear coefficient
lap = r.Dx.^2+r.Dy.^2;                         %%Laplacian
L(1,:) = 1i*0.5*lap(:);                        %%Linear coefficient(1)
L(2,:) = 1i*0.5*lap(:);                        %%Linear coefficient(2)
end                                            %end linear coefficient

