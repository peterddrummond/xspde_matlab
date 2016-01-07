function [e] = GaussianAv7D()                                  %%main function

in.name =        'Diffraction in 7 dimensions';
in.dimension =   7; 
in.points =   [5,11,11,11,11,11,11]; %%dimension: 1-4
in.ranges =   [10,6,6,6,6,6,6];
in.initial =     @(w,r) exp(-0.5*(r.x.^2+r.y.^2+r.z.^2+r.x5.^2+r.x6.^2+r.x7.^2));  
in.da =          @(a,~,~) zeros(size(a));                   %%da 
in.observe{1} =  @(a,r) a.*conj(a);                         %%observable 1
in.observe{2} =  @(a,r) xave(a.*conj(a));                   %%observable 2
in.observe{3} =  @(a,r) xint(a.*conj(a),r);                 %%observable 3
in.observe{4} =  @(a,r) a.*conj(a);                         %%observable 4
in.observe{5} =  @(a,r) xave(a.*conj(a));                   %%observable 5
in.observe{6} =  @(a,r) xint(a.*conj(a),r.dk,r);            %%observable 6
in.linear =      @(D,r) 1i*0.05*(D.x.^2+D.y.^2+D.z.^2+D.x5.^2+D.x6.^2+D.x7.^2);
in.transforms =  {0,0,0,[0,1,1,1,1,1,1],[0,1,1,1,1,1,1],[0,1,1,1,1,1,1]};     
in.images =      2;                                         %%number of images
in.transverse =  2;                                         %%transverse plots
in.olabels =    {'I','<I>','\int I dV','I(k)','<I(k)>','\int I dK'};%%labels 
in.compare{1} =  @(t,~) [1+(t/10).^2].^(-3);              %%comparison 
in.compare{4} =  @(t,~) 1+0.*t;                             %%comparison 
in.pdimension = {3,1,1,1,1,1};
e  = xspde(in);                                             %%simulation
end                                                         %%end of main 
