function  [e]   = Sphere10()
% e = SPHERE10()  simulates 10d spherical diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in.name       = 'Diffusion on a 10-dimensional hypersphere';
in.C          = 1;
in.qc         = eye(10);
in.iterproj   = 3;
in.X0         = [1,0,0,0,0,0,0,0,0,0]';
in.fields     = 10;
in.ranges     = 5;
in.points     = 51;
in.ensembles  = [40, 10, 10];
in.project    = @Quadproj;
in.compare{1} = @(r) 2*(1-exp(-4.5*r.t));
in.compare{2} = @(r) 1+0*r.t;
in.deriv      = @(a, w, r)  w;
in.initial    = @(w, r) r.X0;  
in.observe{1} = @(a, r) sum((r.X0-a).^2,1);
in.observe{2} = @(a, r) (abs(sum(a.*(r.qc*a),1)));
in.diffplot   = {1,0};
in.olabels    = {'\langle R^2\rangle','\langle |f|\rangle'};
in.step       = @MPnproj;
e             = xspde(in);   
end
