function [e] = Gaussian()                      %%name of main function

input.name =       'Gaussian';                 %%name for simulation
input.dimension =  4;                          %%dimension: 1-4 = t,x,y,z
input.ranges =     [2,10,10,10];               %%ranges: t,x,y,z
input.graphs =     [1,1];                      %%x,k graphs
input.Initialise = @Initialise;                %%Initialisation  handle
input.Deriv =      @Deriv;                     %%Derivative  handle
input.Observe =    @Observe;                   %%Derivative  handle
input.Linear =     @Linear;                    %%Derivative  handle
graph.images =     [4,2];                      %%number of graphics images
graph.imagetype =  [1,1];                      %%type of graphics images
graph.transverse = [2,2];                      %%transverse plots
graph.olabels =    {'|a(x)|^2','|a(k)|^2'};    %%labels for observables
graph.compares =   [1,1];                      %%comparisons?
graph.Compare =    @Compare;                   %%Comparison handle

[e,data,lattice] = Stochastic(input);          %%Stochastic program
ec = Graphics(graph,data,lattice);             %%Graphics data
e = max(e,ec);
end                                            %%end of main function

%%XSPDE user functions

function a0 = Initialise(w,in)                 %%Initialises fields
[x,y,z] = meshgrid (in.x{2},in.x{3},in.x{4});  %%lattice vectors
a0 = exp(-0.5*(x.^2+y.^2+z.^2));
a0 = reshape(a0,in.dma);                       %%reshape fields
end                                            %%end initialise fields

function da_dt  =  Deriv(a,~,dt,~,in)          %%Derivatives
da_dt  =  zeros(in.dma);                       %%zero derivative
end                                            %%end local derivatives

function L = Linear(in)                        %%Linear coefficient
[kx,ky,kz] = meshgrid (in.k{2},in.k{3},in.k{4});
L = -1i*0.5*(kx.^2+ky.^2+kz.^2);               %%3 transverse dimensions
L=reshape(L,[1,in.points(2:in.dimension)]);
end                                            %end linear coefficient

function o  =  Observe(a,ft,~,~)               %%Observables
o(1,:) = a(1,:).*conj(a(1,:));                 %%intensity
a = Fourier(a,ft);                             %%Fourier transform
o(2,:) = a(1,:).*conj(a(1,:));                 %%mean in k-space
end                                            %%end observables

function c  =  Compare(t,~)                    %%Comparison functions
c(1,:)= [1+t.^2].^(-3/2);                      %%Used if compare(1)=True
c(2,:)= 1.;                                    %%Used if compare(2)=True
end                                            %%end comparisons 