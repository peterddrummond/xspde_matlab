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

p.name        =  'KE of a gaussian';
p.xlabels     =  {'t', 'z', 'x'};
p.klabels     =  {'\omega', 'k_z', 'k_x'};
p.headers     =  {'','',''};
p.font        =  {10,12,14,16,18,20};
p.dimensions  =  3;
p.fields      =  1;
p.ranges      =  [5 200 6];
y0            =  0.235;
x0            =  3.3;
p.points      =  [49 70 34];
p.initial     =  @(v,p) sqrt(exp(-p.x.^2/x0^2 - p.y.^2/y0^2)/(pi*x0*y0));
p.images      =  {2 2};
p.pdimension  =  {3 3 1 1 1 1 1 1};
p.linear      =  @(p) 1i*(p.Dx.^2);
p.observe{1}  =  @(a,p)  abs(a).^2;
p.transforms  =  {0,[0 1 1],[0 1 1],0,0,[0 1 1]};
p.observe{2}  =  @(a,p)  abs(a).^2;
p.observe{3}  =  @(a,p)  Int(abs(a).^2.*(p.kx.^2 + p.ky.^2), p.dk, p);
p.compare{3}  =  @(p) 1/(2*x0^2) + 1/(2*y0^2);
p.observe{4}  =  @(a,p) sqrt(Int(p.y.^2.*abs(a).^2,p.dx,p)./Int(abs(a).^2,p.dx,p));
p.compare{4}  =  @(p) y0/sqrt(2);
p.observe{5}  =  @(a,p) sqrt(Int(p.x.^2 .*abs(a).^2,p.dx,p)./Int(abs(a).^2,p.dx,p));
p.observe{6}  =  @(a,p) sqrt(Int(p.ky.^2.*abs(a).^2,p.dk,p)./Int(abs(a).^2,p.dk,p));
p.compare{6}  =  @(p) 1/(sqrt(2)*y0);
p.olabels     =  {'<|\psi|^2>', '<|\psi(k)|^2>', '<T>', '\sigma_x','\sigma_z','\sigma_{kx}'} ;
e             =  xspde(p);
end	% function test
