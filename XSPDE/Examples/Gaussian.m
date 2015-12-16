function [e] = Gaussian()                                  %%main function

in.name =       'Gaussian diffraction';
in.dimension =  4;                                         %%dimension: 1-4
in.initial =    @(w,r) exp(-0.5*(r.x.^2+r.y.^2+r.z.^2));   %%initialisation  
in.da =         @(a,~,~) zeros(size(a));                   %%da 
in.observe{1} = @(a,~) a.*conj(a);                         %%x-observable
in.observe{2} = @(a,~) a.*conj(a);                         %%k-observable
in.linear =     @(D,r) 1i*0.05*(D.x.^2+D.y.^2+D.z.^2);     %%laplacian
in.transforms = {0,[0,1,1,1]};
in.images =     {2};                                       %%number of images
in.transverse = {1};                                       %%transverse plots
in.olabels =    {'|a(r)|^2','|a(k)|^2'};                   %%labels
in.compare =    {@(t,~) [1+(t/10).^2].^(-3/2)};            %%comparison fn 
e  = xspde(in);                                            %%simulation
end                                                        %%end of main 
