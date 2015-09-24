function [e] = Gain()                          %%name of main function

%%XSPDE input parameters

input.name =       'Loss with noise';          %%name for simulation
input.ranges =     4;                          %%ranges: t,x,y,z  
input.noises =     [2,0];                      %%xnoises, knoises per point
input.ensembles =  [100,16,1];                  %%samples,ensembles,parallel
input.Initialise = @Initialise;                %%Initialisation  handle
input.Deriv =      @Deriv;                     %%Derivative  handle
input.Observe =    @Observe;                   %%Observe  handle
graph.headers =    [1,1];                      %%headers in graphs
graph.olabels =    {'|a|^2'};                  %%labels for observables
graph.compares =   1;                          %%comparisons?
graph.Compare =    @Compare;                   %%Comparison handle
input2        =    input;                      %%Second inputs
input2.steps  =     4;                         %%Steps per plotted point
input2.name =       'Gain with noise';         %%name for simulation
input2.Deriv =      @Deriv2;                   %%Derivative  handle
graph2 = graph;                                %%Second graphs
graph2.Compare =   @Compare2;                  %%Comparison handle


[e,data,lattice] = Stochastic({input,input2}); %%Stochastic program
ec =Graphics({graph,graph2},data,lattice);     %%Graphics program
e = max(e,ec);                                 %%Max errors
end                                             
%%end of main function

%%XSPDE user functions

function a0 = Initialise(~,in)                 %%Initialises fields
a0 = ones(in.dma);                             %%initial value
end                                            %%end initialise fields

function da  =  Deriv2(a,~,dt,dw,~)            %%Derivatives
da(1,:)  = a(1,:)*dt + dw.x(1,:)+1i*dw.x(2,:); %%Gain equation
end                                            %%end local derivatives

function da  =  Deriv(a,~,dt,dw,~)             %%Derivatives
da(1,:) = -a(1,:)*dt + dw.x(1,:)+1i*dw.x(2,:); %%Gain equation
end                                            %%end local derivatives

function o  =  Observe(a,~,~,~)                %%Observables
o(1,:) = a(1,:).*conj(a(1,:));                 %%field
end                                            %%end observables

function c  =  Compare2(t,~)                   %%Comparison functions
c = 2*exp(2*t)-1;                              %%used if compare(1)=True
end                                            %%end comparisons 

function c  =  Compare(t,~)                    %%Comparison functions
c = 1+0*t;                                     %%used if compare(1)=True
end                                            %%end comparisons 