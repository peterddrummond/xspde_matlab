function [e] = SHO()                           %%name of main function

%%XSPDE simulation code

input.ranges =     2*pi;                       %%ranges: t,x,y,z  
input.Initialise = @Initialise;                %%Initialisation  handle
input.Deriv =      @Deriv;                     %%Derivative  handle
graph.olabels =    {'\Re(a)'};                 %%labels for observables
graph.compares =   1;                          %%comparisons?
graph.Compare =      @Compare;                 %%Derivative  handle

[e,data,lattice] = Stochastic(input);          %%Stochastic program
ec =Graphics(graph,data,lattice);              %%Graphics program
e =max(e,ec);                                  %%errors
end                                            %%end of main function

%%XSPDE user functions

function a0 = Initialise(~,in)                 %%Initialises fields
a0 = ones(in.dma);                             %%initial value
end                                            %%end initialise fields

function da  =  Deriv(a,~,dt,~,~)              %%Derivatives
da  =  1i*a*dt;                                %%SHO equation
end                                            %%end local derivatives

function c  =  Compare(t,in)                   %%Comparison functions
c = cos(t);                                    %%used if compare(1)=True
end                                            %%end comparisons 