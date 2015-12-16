function [e] = GaussianAv1D()                             %%main function

in.name =        'Gaussian diffraction in 1D';
in.dimension =   2;                                       %%dimension: 1-4
in.initial =     @(w,r) exp(-0.5*(r.x.^2));               %%initialisation  
in.da =          @(a,~,~) zeros(size(a));                 %%da 
in.observe{1} =  @(a,r) a.*conj(a);                       %%observable 1
in.observe{2} =  @(a,r) xave(a.*conj(a),r.dx,r);          %%observable 2
in.observe{3} =  @(a,r) xint(a.*conj(a),r.dx,r);          %%observable 3
in.observe{4} =  @(a,r) a.*conj(a);                       %%observable
in.linear =      @(D,in) 1i*0.1*(D.x.^2);                 %%laplacian
in.transforms =  {0,0,0,[0,1]};                           %%transforms
in.transverse =  {1};                                     %%transverse plots
in.olabels =     {'I','<I>','\int I dV','I(k)'};          %%labels 
in.compare{1} =  @(t,~,~) [1+(t/5).^2].^(-1/2);           %%comparison 
in.compare{4} =  @(t,~,~) 1+0.*t;                         %%comparison 
e  = xspde(in);                                           %%simulation
end                                                       %%end of main 
