function [e] = Gaussian()                                  %%main function

in.name =       'Gaussian diffraction with loss';
in.dimension =  4;                                         %%dimension: 1-4
in.initial =    @(w,r) exp(-0.5*(r.x.^2+r.y.^2+r.z.^2));   %%initialisation  
in.da =         @(a,~,~) zeros(size(a));                   %%da 
in.observe{1} = @(a,~) a.*conj(a);                         %%x-observable
in.linear =     @(D,r) -0.5+1i*0.05*(D.x.^2+D.y.^2+D.z.^2);%%laplacian
in.transforms = {0};
in.images =     {2};                                       %%number of images
in.transverse = {1};                                       %%transverse plots
in.olabels =    {'|a(r)|^2'};                   %%labels
in.compare =    @(t,~) [1+(t/10).^2].^(-3/2).*exp(-t);    %%comparison fn 
e  = xspde(in);                                            %%simulation
end                                                        %%end of main 
