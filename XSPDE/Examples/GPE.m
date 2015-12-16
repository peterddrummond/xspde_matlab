function [e] = GPE()                                    %%name of function
%% Simulates a noisy GPE in a harmonic well with an initial vortex
cd ~
in.name =  'GPE';
in.dimension =  3;                                     
in.points = [201,128,128];
in.ranges = [1,20,20];
in.noises = [0,2];
in.rfilter =    @(r)          [exp(- r.kx.^2-r.ky.^2);exp(- r.kx.^2-r.ky.^2)];
in.nfilter =    @(r)          [exp(- r.kx.^2-r.ky.^2);exp(- r.kx.^2-r.ky.^2)];
in.initial =    @(w,r)        (r.x+i*r.y)./(1+10*(r.x.^2 +r.y.^2));
V =             @(r)          0.01*(r.x.^2 +r.y.^2)-0.001*i*(r.x.^2 +r.y.^2).^2;
b =             @(xi)         .1*(xi(1,:)+i*xi(2,:));
in.da =         @(a,xi,r)     -i*a.*(V(r)+conj(a).*a)+b(xi);
in.linear =     @(D,r)         0.5*i*(D.x.^2+D.y.^2);
in.observe{1} = @(a,r)         a.*conj(a);
in.images  =  {2};                                       
in.imagetype   = {2}; 
in.olabels = {'|a|^2'};
in.file =       'GPE.h5';
[e,in] = xsim(in);                         %%sim program
xgraph(in);                                %%graph program
e = max(e);
end                                        %%end of function
