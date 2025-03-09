function [error] = Kubo()

%%XSPDE parameters

input.name =       'Kubo Oscillator';          %%name for simulation
input.ranges =     5;                          %%ranges: t,x,y,z 
input.fields =     2;                          %%field number 
input.ensembles =  [400,16,1];                 %%samples,ensembles,parallel
input.noises =     [2,0];                      %%noises in x and k
input.Initialise = @Initialise;                %%Initialisation  handle
input.Deriv =      @Deriv;                     %%Derivative  handle
input.Matlabfile = 'Kubo';                     %%Test file output
graph.olabels =    {'<a_1>','<a_2>'};          %%labels for observables
graph.compares =   [1,1];                      %%comparisons?
graph.Compare =    @Compare;                   %%Comparison handle

[error,data,input] = Stochastic(input);        %%Stochasic program
differences = Graphics(graph,data,input);      %%Graphics data
error =max(error,differences);                 %%Error max
end

function a0 = Initialise(~,in)                 %%Initialises fields
a0 = ones(in.dma);                             %%initial value
end                                            %%end initialise fields

function da  =  Deriv(a,~,dt,dw,in)            %%Derivatives
da(1,:)  =  1i*dw.x(1,:).*a(1,:);              %%Kubo equation
da(2,:)  =  1i*dw.x(2,:).*a(2,:);              %%Kubo equation 2
da(2,:)  =  da(2,:) +  1i*dt*a(2,:);           %%Kubo equation + detuning
end                                            %%end local derivatives

function c  =  Compare(t,~)                    %%Comparison functions
c(1,:) = exp(-t/2);                            %%used if compare(1)=True
c(2,:) = cos(t).*exp(-t/2);                    %%used if compare(2)=True
end                                            %%end comparisons 