function [e] = GaussianKinetic()
%   e  =  GAUSSIANKINETIC() tests xSPDE for a linear Schroedinger equation.
% This computes the kinetic energy of a steady Gaussian wavefunction, as a 
% test of the code for calculating moments in position and momentum space. 
% It also demonstrates how to change the  graphics labels on axes, 
% use of variables inside a function,  and images and pdimension inputs.
% Axes are deliberately labeled differently from the internal labels,
% and the code calculates the uncertainty in position and momentum
%   Tests a (two+one)-dimensional partial differential equation for:
%   (1) Initial asymmetric gaussian in space
%   (2) Using constants defined in the top function definition
%   (3) Using transforms  for some, not all observables
%   (4) Integrating over powers of momentum in the momentum grid
%   (5) Suppressing headers in first three graphs, and changing fonts
%   (6) Using transforms over part of the spatial grid
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

in.name =	    'KE of a gaussian';
in.xlabels =    {'t', 'z', 'x'};
in.klabels =    {'\omega', 'k_z', 'k_x'};
in.headers =   {'','',''};
in.font =      {10,12,14,16,18,20};
in.dimension =	3;
in.fields =		1; 
in.ranges =		[5 200 6];	
y0    = 		0.235;
x0     = 		3.3;
in.points =		[49 70 34];
in.initial =    @(v,r) sqrt(exp(-r.x.^2/x0^2 - r.y.^2/y0^2)/(pi*x0*y0));
in.images =		{2 2};
in.pdimension =	{3 3 1 1 1 1 1 1};
in.linear =      @(r) 1i*(r.Dx.^2);
in.observe{1} =	@(a,r)  abs(a).^2;
in.transforms = {0,[0 1 1],[0 1 1],0,0,[0 1 1]};
in.observe{2} =	@(a,r)  abs(a).^2;
in.observe{3} =	@(a,r)  xint(abs(a).^2.*(r.kx.^2 + r.ky.^2), r.dk, r);
in.compare{3} = @(t, in) 1/(2*x0^2) + 1/(2*y0^2);
in.observe{4} = @(a,r) sqrt(xint(r.y.^2.*abs(a).^2,r.dx,r)./xint(abs(a).^2,r.dx,r));
in.compare{4} = @(t,in) y0/sqrt(2);
in.observe{5} = @(a,r) sqrt(xint(r.x.^2.*abs(a).^2,r.dx,r)./xint(abs(a).^2,r.dx,r));
in.observe{6} = @(a,r) sqrt(xint(r.ky.^2.*abs(a).^2,r.dk,r)./xint(abs(a).^2,r.dk,r));
in.compare{6} = @(t,in) 1/(sqrt(2)*y0);
in.olabels =	{'<|\psi|^2>', '<|\psi(k)|^2>', '<T>', '\Delta x','\Delta z','\Delta k_x'} ;
e             = xspde(in);
end	% function test
