function [e] = SolitonAv()                                

in.name =          'NLS soliton with integration';
in.dimension =     2;                                     %%dimension: 1-4  
in.print =         0;
in.observe{1} =    @(a,r)      xave(a,r.dx,r);
in.observe{2} =    @(a,r)      xint(a,r.dx,r);
in.initial =       @(w,r)      sech(r.x);                 %%Initialisation 
in.da =            @(a,~,r)    i*a.*(conj(a).*a);         %%Derivative 
in.linear =        @(D,r)      0.5*i*(D.x.^2-1.0);        %%Laplacian
in.compare   =     {@(t,~) pi/10+0*t,@(t,~) pi+0*t};
in.olabels   =     {'<<a(x)>>','\int<a(x)>dx'};  %%Labels
in.pdimension =    [1,1];
[e,data] =         xspde(in);                          %%main program
end                                                       %%end of function