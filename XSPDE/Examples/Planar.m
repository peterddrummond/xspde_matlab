function [e] = Planar()                        %%name of main function

in.name =       'Planar noise growth';         %%name for simulation
in.dimension =  3;                             %%dimension: 1-4 = t,x,y,z
in.fields =     2;                             %%field components
in.ranges =     [1,5,5];                       %%ranges: t,x,y,z  
in.steps =      2;                             %%steps per plotted point
in.noises =     [2,2];                         %%xnoises, knoises per point
in.ensembles =  [10,2,2];                      %%samples,ensembles,parallel
in.initial =    @Initial;                      %%Initialisation  handle
in.da  =        @Da;                           %%Derivative  handle
in.linear =     @Linear;                       %%Derivative  handle

in.observe{1} = @(a,~) xave(a(1,:).*conj(a(1,:)));%%Observe  handle
in.observe{2} = @(a,~) xave(a(2,:).*conj(a(2,:)));%%Observe  handle
in.observe{3} = @(a,~) xave(a(1,:).*conj(a(2,:)));%%Observe  handle
in.transforms = {[0,0,0],[0,1,1],[0,1,1]};

in.olabels{1} = '<<|a_1(x)|^2>>';              %%labels  
in.olabels{2} = '<<|a_2(k)|^2>>';              %%labels
in.olabels{3} = '<<a_1(k)a^*_2(k)>>';          %%labels 
in.compare{1} = @(t,in) [1+t]/in.dV;
in.compare{2} = @(t,in) [1+t]/in.dK;
in.compare{3} = @(t,in) 0*t;
in.pdimension = {1,1,1,1};                     %%maximum plot dimension

e  =  xspde(in);                            %%Stochasic program
end                                            %%end of main function

%%XSPDE user functions

function a0 = Initial(w,r)                  %%Initialises fields
a0(1,:) = (w(1,:)+1i*w(2,:))/sqrt(2);          %%lattice vectors
a0(2,:) = (w(3,:)+1i*w(4,:))/sqrt(2);          %%lattice vectors
end                                            %%end initialise fields

function da  =  Da(a,xi,r)                  %%Derivatives

da(1,:) = (xi(1,:)+1i*xi(2,:))/sqrt(2);        %%complex k-noise equation
da(2,:) = (xi(3,:)+1i*xi(4,:))/sqrt(2);        %%complex x-noise equation
end                                            %%end local derivatives

function L = Linear(D,r)                      %%Linear coefficient
lap = D.x.^2+D.y.^2;                           %%Laplacian
L(1,:) = 1i*0.5*lap(:);                        %%Linear coefficient(1)
L(2,:) = 1i*0.5*lap(:);                        %%Linear coefficient(2)
end                                            %end linear coefficient

