function  [e]   = Sphere10()
% e = SPHERE10()  simulates 10d spherical diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.name       = 'Sphere10: Diffusion on a 10D hypersphere';
p.C          = 1;
p.qc         = eye(10);
p.iterproj   = 3;
p.X0         = [1,0,0,0,0,0,0,0,0,0]';
p.fields     = 10;
p.ranges     = 5;
p.points     = 51;
p.ensembles  = [40, 10, 10];
p.project    = @Quadproj;
p.compare{1} = @(r) 2*(1-exp(-4.5*r.t));
p.compare{2} = @(r) 1+0*r.t;
p.deriv      = @(a, w, r)  w;
p.initial    = @(w, r) r.X0;  
p.observe{1} = @(a, r) sum((r.X0-a).^2,1);
p.observe{2} = @(a, r) (abs(sum(a.*(r.qc*a),1)));
p.diffplot   = {1,0};
p.olabels    = {'< R^2 >','< |f| >'};
p.method     = @MPnproj;
e            = xspde(p);   
end
