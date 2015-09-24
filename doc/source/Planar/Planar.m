function [e] = Planar()                        %%name of main function

%%XSPDE input parameters

input.name =       'Planar';                   %%name for simulation
input.dimension =  3;                          %%dimension: 1-4 = t,x,y,z
input.fields =     2;                          %%field components
input.ranges =     [1,5,5];                    %%ranges: t,x,y,z  
input.steps =      2;                          %%steps per plotted point
input.noises =     [2,2];                      %%xnoises, knoises per point
input.ensembles =  [10,16,1];                  %%samples,ensembles,parallel
input.graphs =     [2,2];                      %%graph functions in x and k
input.Initialise = @Initialise;                %%Initialisation  handle
input.Deriv =      @Deriv;                     %%Derivative  handle
input.Linear =     @Linear;                    %%Derivative  handle
input.Observe =    @Observe;                   %%Observe  handle

[e,data,lattice] = Stochastic(input);          %%Stochastic program

%%XSPDE graphics parameters

graph.images =     [4,2,0,0];                  %%number of graphics images
graph.imagetype =  [1,1,1,1];                  %%type of graphics images
graph.transverse = [2,2,0,0];                  %%transverse plots
graph.headers =    [1,1,1,1];                  %%headers in graphs
graph.pdimension = [4,1,1,1];                  %%maximum plot dimension
graph.font =       18;                         %%label fontsize
graph.olabels =    {'<|a_1(x)|^2>','<<|a_1(x)|^2>>',...
    '<<|a_2(k)|^2>>','<<a_1(k)a^*_2(k)>>'};    %%labels
graph.compares =   [1,1,1,1];                  %%comparisons?
graph.Compare =    @Compare;                   %%Comparison handle

ec = Graphics(graph,data,lattice);             %%Graphics data
e = max(ec,e);
end                                            %%end of main function

%%XSPDE user functions

function a0 = Initialise(w,in)                 %%Initialises fields
w1 = Nifft(w.k);                               %%inverse FFT of k-noise
w = w.x;                                       %%x-noise
a0(1,:) = (w1(1,:)+1i*w1(2,:))/sqrt(2);        %%lattice vectors
a0(2,:) = (w(1,:)+1i*w(2,:))/sqrt(2);          %%lattice vectors
end                                            %%end initialise fields

function da  =  Deriv(a,t,dt,dw,in)            %%Derivatives
dw1 = Nifft(dw.k);                             %%inverse FFT of k-noise
dw =dw.x;                                      %%x-noise
da(1,:) = (dw1(1,:)+1i*dw1(2,:))/sqrt(2);      %%complex k-noise equation
da(2,:) = (dw(1,:)+1i*dw(2,:))/sqrt(2);        %%complex x-noise equation
end                                            %%end local derivatives

function L = Linear(in)                        %%Linear coefficient
[kx,ky] = meshgrid (in.k{2},in.k{3});          %%initialize k on grid
L = zeros(in.dml);                             %%initialize L
lap = -(kx.^2+ky.^2);                          %%Laplacian
L(1,:) = 1i*0.5*lap(:);                        %%Linear coefficient(1)
L(2,:) = 1i*0.5*lap(:);                        %%Linear coefficient(2)
end                                            %end linear coefficient

function o  =  Observe(a,ft,~,~)               %%Observables
o(1,:) = a(1,:).*conj(a(1,:));                 %%intensity
o(2,:) = mean (o(1,:));                        %%mean intensity
a = Fourier(a,ft);                             %%Fourier transform
o(3,:) = mean (a(2,:).*conj(a(2,:)));          %%mean in k-space
o(4,:) = real(mean (a(1,:).*conj(a(2,:))));    %%correlation in k-space
end                                            %%end observables

function c  =  Compare(t,in)                   %%Comparison functions
V = in.ranges(2)*in.ranges(3);
N = in.points(2)*in.points(3);
c(1,:)= [1+t]*(N/V);                           %%Used if compare(1)=True
c(2,:)= c(1,:);                                %%Used if compare(2)=True
c(3,:)= [1+t]*V/(2*pi)^2;                      %%Used if compare(3)=True
c(4,:)= 0;                                     %%Used if compare(4)=True
end                                            %%end comparisons 