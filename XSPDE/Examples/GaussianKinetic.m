function [e] = GaussianKinetic()                        
% This computes the kinetic energy of a steady Gaussian wavefunction, as a 
% test of the code for calculating moments in position and momentum space. 
% It also demonstrates how to change the  graphics labels on axes, 
% use of variables inside a function,  and images and pdimension inputs.
% Axes are deliberately labeled differently from the internal labels,
% and the code calculates the uncertainty in position and momentum

in.xlabels =    {'t', 'z', 'x'};
in.klabels =    {'\omega', 'k_z', 'k_x'};
tmax =          5;
in.name =	    'KE of a gaussian';
in.dimension =	3;
in.fields =		1; 
in.ranges =		[tmax 200 6];	
y0    = 		0.235;
x0     = 		3.3;
in.points =		[49 70 34];
in.initial =    @(v,r) sqrt(exp(-r.x.^2/x0^2 - r.y.^2/y0^2)/(pi*x0*y0));
in.da =		    @(a,~,r) zeros(size(r.x));
in.ensembles =	[1 1 1];
in.images =		{2 2 0 0 0 0 0 0};
in.pdimension =	{3 3 1 1 1 1 1 1};
in.observe{1} =	@(a,r)  abs(a).^2;
in.olabels{1} =	'<|\psi|^2>';
in.transforms{2} = [false true true];
in.observe{2} =	@(a,r)  abs(a).^2;
in.olabels{2} =	'<|\psi(k)|^2>';
in.transforms{3} = [false true true];
in.observe{3} =	@(a,r)  xint(abs(a).^2.*(r.kx.^2 + r.ky.^2), r.dk, r);
in.compare{3} = @(t, in) 1/(2*x0^2) + 1/(2*y0^2);
in.olabels{3} =	'<T>';
in.olabels{4} = '\Delta x';
in.observe{4} = @(a,r) sqrt(xint(r.y.^2.*abs(a).^2,r.dx,r)./xint(abs(a).^2,r.dx,r));
in.compare{4} = @(t,in) y0/sqrt(2);
in.olabels{5} = '\Delta z';
in.observe{5} = @(a,r) sqrt(xint(r.x.^2.*abs(a).^2,r.dx,r)./xint(abs(a).^2,r.dx,r));
in.compare{5} = @(t,in) x0/sqrt(2);
in.transforms{6} = [false true true];
in.olabels{6} = '\Delta k_x';
in.observe{6} = @(a,r) sqrt(xint(r.ky.^2.*abs(a).^2,r.dk,r)./xint(abs(a).^2,r.dk,r));
in.compare{6} = @(t,in) 1/(sqrt(2)*y0);
e             = xspde(in);
end	% function test
