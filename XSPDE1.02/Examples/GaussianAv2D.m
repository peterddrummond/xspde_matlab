function [e] = GaussianAv2D()                             %%main function

in.name =        'Gaussian diffraction in 2D';
in.dimension =   3;                                       %%dimension: 1-4
in.initial =     @(w,r) exp(-0.5*(r.x.^2+r.y.^2));        %%initialisation  
in.da =          @(a,~,~) zeros(size(a));                 %%da 
in.observe{1} =  @(a,r) a.*conj(a);                       %%observable 1
in.observe{2} =  @(a,r) xave(a.*conj(a),r.dx,r);          %%observable 2
in.observe{3} =  @(a,r) xint(a.*conj(a),r.dx,r);          %%observable 3
in.observe{4} =  @(a,r) a.*conj(a);                       %%observable 4
in.linear =      @(D,r) 1i*0.05*(D.x.^2+D.y.^2);          %%laplacian
in.transforms =  {0,0,0,[0,1,1]};                         %%transforms
in.images =      {2};                                     %%number of images
in.transverse =  {2};                                     %%transverse plots
in.olabels =     {'I','<I>','\int I dV','I(k)','<I(k)>','\int I dK'};%%labels 
in.compare{1} =  @(t,~,~) [1+(t/10).^2].^(-1);            %%comparison 
in.pdimension =  {3,1,1,1};
e  = xspde(in);                                           %%simulation
end                                                       %%end of main 
