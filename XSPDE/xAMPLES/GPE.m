function [e] = GPE()                                   
%   e  =  GPE() tests xSPDE for  a vortex in two space dimensions
%   Tests a three-dimensional noisy GPE in a harmonic well by
%   (1) Modifying points and ranges in 3 dimensions
%   (2) Using momentum filters on noises generated
%   (3) Setting the noise and random transform filter in k-space
%   (4) Changing the movie image type to give a contour plot
%   (5) Using auxiliary functions defined inline
%   (6) Graphing data from results saved in an HDF5 file
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

cd ~
in.name =  'GPE';
in.dimension =  3;                                     
in.points = [101,64,64];
in.ranges = [1,20,20];
in.noises = [0,2];
in.rfilter =    @(r)          [exp(- r.kx.^2-r.ky.^2);exp(- r.kx.^2-r.ky.^2)];
in.nfilter =    @(r)          [exp(- r.kx.^2-r.ky.^2);exp(- r.kx.^2-r.ky.^2)];
b =             @(xi)         .1*(xi(1,:)+1i*xi(2,:));
in.initial =    @(w,r)        (r.x+1i*r.y)./(1+10*(r.x.^2 +r.y.^2))+b(w);
V =             @(r)          0.01*(r.x.^2 +r.y.^2)-0.001*1i*(r.x.^2 +r.y.^2).^2;
in.da =         @(a,xi,r)     -1i*a.*(V(r)+conj(a).*a)+b(xi);
in.linear =     @(r)          0.5*1i*(r.Dx.^2+r.Dy.^2);
in.observe{1} = @(a,r)        a.*conj(a);
in.images  =  {2};                                       
in.imagetype   = {2}; 
in.olabels = {'|a|^2'};
in.file =       'GPE.mat';
in.print =2;
[e,in] = xsim(in);                         %%sim program
xgraph(in.file,in);                        %%graph program
end                                        %%end of function