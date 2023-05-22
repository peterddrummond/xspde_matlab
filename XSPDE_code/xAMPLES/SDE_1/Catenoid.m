function  [e]       =  Catenoid
%   e  =  CATENOID() tests xSPDE for a catenoid diffusion
%   Tests for:
%   (1) Catenoid projections
%   (2) Difference plots
%   (3) Using a projected algorithm
%   (4) Omitting an observe function
%   (5) Omitting a compare function
%   (6) Omitting a diffplot
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License

p.iterproj      =  3;
p.X0            =  [1,0,0]';
p.fields        =  3;
p.ranges        =  5;
p.points        =  51;
p.steps         =  1;
p.ensembles     =  [400, 10];
p.compare{2}    =  @(p) 2*p.t;
p.deriv         =  @(a, w, p)  w;
p.initial       =  @(w, p) p.X0 + 0.*w;
p.observe{2}    =  @(a, p) sum((p.X0-a).^2,1);
p.diffplot{2}   =  1;
p.function{1}   =  @(o, p) o{2}.^2;
p.olabels       =  {'< R^2 >','< R^2 >'};
p.project       =  @Catproj;
p.name          =  'Catenoid diffusion in 3D';
p.method        =  @MPnproj;
e               =  xspde(p);
end
