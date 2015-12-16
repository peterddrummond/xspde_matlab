function [e] = GPE2()                                    %%name of function
%% Simulates a 1D GPE in a harmonic well with two fields
in.name =  'GPE2';
in.dimension =  2; 
in.fields =  2;  
in.points = [101,64];
in.ranges = [20,20];
in.steps = 4;
in.noises = [2,0];
in.initial =    @(w,r)          [sech(r.x+5);sech(r.x-5)];
V        =      @(r)            0.1*r.x.^2;
rho        =    @(a)            conj(a(1,:)).*a(1,:)+conj(a(2,:)).*a(2,:);
da1        =    @(a,r)           -i*a(1,:).*(V(r)+rho(a));
da2        =    @(a,r)           -i*a(2,:).*(V(r)+rho(a));
in.da =         @(a,xi,r)       [da1(a,r);da2(a,r)];
in.linear =     @(D,r)          [0.5*i*(D.x.^2+D.y.^2);0.5*i*(D.x.^2+D.y.^2)]; 
in.observe{1} = @(a,r) a(1,:).*conj(a(1,:));
in.observe{2} = @(a,r) a(2,:).*conj(a(2,:));
in.images  =  {2};                                       
in.imagetype   = {2};  
in.olabels = {'|a_1|^2','|a_2|^2'};
[e,data] = xspde(in);                             %%main program
end                                                     %%end of function
