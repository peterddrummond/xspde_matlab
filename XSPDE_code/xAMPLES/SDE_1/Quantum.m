function e      =  Quantum()
%   e  =  QUANTUM() tests xSPDE for an SDE with a known spectrum
%   This is the quantum harmonic oscillator with a vacuum input and
%   output, for which one expects 0.5 symmetrically ordered 'photons'
%   per mode.

p.name          =  'Quantum spectrum';               %%name for simulation
p.points        =  601;                              %%points in time
p.ranges        =  160;                              %%range in time
p.order         =  2;
p.auxfields     =  2;
p.noises        =  2;                                %%noises per point
p.ensembles     =  [200,1,12];                       %%samples,ensembles
p.initial       =  @(w,~) (w(1,:)+1i*w(2,:))/(2);    %%Initialisation
p.a1            =  @(w) (w(1,:)+1i*w(2,:))/2;
p.deriv         =  @(a,w,~)  -a(1,:)+sqrt(2)*p.a1(w);%%Derivative
p.define        =  @(a,w,p) [p.a1(w);sqrt(2)*a(1,:)-p.a1(w)];
T               =  @(p) p.ranges(1);
p.observe{1}    =  @(a,p) (2.*pi/T(p))*a(1,:).*conj(a(1,:));
p.observe{2}    =  @(a,p) (2.*pi/T(p))*a(2,:).*conj(a(2,:));
p.observe{3}    =  @(a,p) (2.*pi/T(p))*a(3,:).*conj(a(3,:));
p.observe{4}    =  @(a,p) a(1,:).*conj(a(1,:));
p.transforms    =  {1,1,1};                          %%Fourier transforms
p.olabels{1}    =  '|a(\omega)|^2';
p.olabels{2}    =  '|a_{in}(\omega)|^2';
p.olabels{3}    =  '|a_{out}(\omega)|^2';
p.olabels{4}    =  '|a|^2';
p.compare{1}    =  @(p) 1./(1+p.w.^2);               %%Comparisons
p.compare{2}    =  @(p) 0.5;                         %%Comparisons
p.compare{3}    =  @(p) 0.5;                         %%Comparisons
p.compare{4}    =  @(p) 0.5;                         %%Comparisons
e               =  xspde(p);                         %%Stochastic toolbox
end
