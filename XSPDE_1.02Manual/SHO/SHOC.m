function [e,data,lattice] = SHO()              %%name of main function

%%XSPDE simulation code

in.ranges =     2*pi;                          %%ranges: t,x,y,z  
in.Initialise = @Initialise;                   %%Initialisation  handle
in.Deriv =      @Deriv;                        %%Derivative  handle
g.olabels =    {'\Re(a)'};                     %%labels for observables
g.compares =   1;                              %%comparisons?
g.Compare =      @Compare;                     %%Derivative  handle

[e,data,lattice] = Stochastic({in});           %%Stochastic program
Graphics({g},data,lattice);                    %%Graphics program
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