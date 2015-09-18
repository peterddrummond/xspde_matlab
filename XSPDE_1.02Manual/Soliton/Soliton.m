function [e] = Soliton()                      %%name of main function

input.name =       'Soliton';                  %%name for simulation
input.dimension =  2;                          %%dimension: 1-4 = t,x,y,z
input.ranges =     [2,20];                     %%ranges: t,x,y,z  
input.Initialise = @Initialise;                %%Initialisation  handle
input.Deriv =      @Deriv;                     %%Derivative  handle
input.Linear =     @Linear;                    %%Derivative  handle
graph.olabels =    {'\Re(a)'};                 %%labels for observables
graph.compares =   1;                          %%comparisons?
graph.Compare =    @Compare;                   %%Comparison handle

[e,data,lattice] = Stochastic(input);          %%Stochastic program
ec = Graphics(graph,data,lattice);             %%Graphics data
e = max(e,ec);                                 %%Combined errors
end                                            %%end of main function

function a0 = Initialise(w,in)                 %%Initialises fields
x = in.x{2};                                   %%lattice vectors
a0 = sech(x);
a0 = reshape(a0,in.dma);                       %%reshape fields
end                                            %%end initialise fields
function da  =  Deriv(a,~,dt,~,~)              %%Derivatives
da  =  1i*a.*(conj(a).*a)*dt;                  %%NLS equation
end                                            %%end local derivatives
function L = Linear(in)                        %%Linear coefficient
kx=in.k{2};                                    %%k-vector
L(1,:) = -.5*1i*(1.0+kx.^2);                   %%laplacian in k-space
end                                            %end linear coefficient
function c  =  Compare(t,~)                    %%Comparison functions
c(1,:)= 1;                                     %%Used if compare(1)=True
end                                            %%end comparisons 