function [e] = SolitonAv()
%   e  =  SOLITONAV() tests xSPDE for a nonlinear Schroedinger equation
%   Tests a two-dimensional partial differential equation for:
%   (1) Putting print to 0 to minimize output printing
%   (2) Averaging over the spatial lattice
%   (3) Setting up the linear response for the interaction picture
%   (4) Compare functions for a 2D simulation
%   (5) Graphing results in 1D and 2D
%   (6) Default values for transverse range and points in x
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =          'NLS soliton with integration';
in.dimension =     2;                                     %%dimension: 1-4  
in.print =         0;
in.observe{1} =    @(a,r)      xave(a,r.dx,r);
in.observe{2} =    @(a,r)      xint(a,r.dx,r);
in.initial =       @(w,r)      sech(r.x);                 %%Initialisation 
in.da =            @(a,~,r)    i*a.*(conj(a).*a);         %%Derivative 
in.linear =        @(r)      0.5*i*(r.Dx.^2-1.0);        %%Laplacian
in.compare   =     {@(t,~) pi/10+0*t,@(t,~) pi+0*t};
in.olabels   =     {'<<a(x)>>','\int<a(x)>dx'};  %%Labels
in.pdimension =    {1,1};
[e,data] =         xspde(in);                          %%main program
end                                                       %%end of function