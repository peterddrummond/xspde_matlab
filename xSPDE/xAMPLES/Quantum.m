function e = Quantum()
%   e  =  QUANTUM() tests xSPDE for an SDE with a known spectrum
%   This is the quantum harmonic oscillator with a vacuum input and 
%   output, for which one expects 0.5 symmetrically ordered 'photons'
%   per mode.

in.name =      'Quantum spectrum';                %%name for simulation
in.points =    601;                               %%points in time
in.ranges =    160;                                %%range in time
in.fields  =   [1,2];                            
in.noises =    2;                                 %%xnoises per point
in.ensembles = [400,1,6];                         %%samples,ensembles
in.initial =   @(w,~) (w(1,:)+1i*w(2,:))/(2);     %%Initialisation
in.a1 =        @(z) (z(1,:)+1i*z(2,:))/2;
in.da =        @(a,z,~)  -a(1,:)+sqrt(2)*a(2,:);  %%Derivative
in.define =    @(a,z,r) [r.a1(z);sqrt(2)*a(1,:)-r.a1(z)];  
T          =   @(r) r.ranges(1);      
in.observe{1}= @(a,r) (2.*pi/T(r))*a(1,:).*conj(a(1,:));             
in.observe{2}= @(a,r) (2.*pi/T(r))*a(2,:).*conj(a(2,:));            
in.observe{3}= @(a,r) (2.*pi/T(r))*a(3,:).*conj(a(3,:)); 
in.observe{4}= @(a,r) a(1,:) - sqrt(2)*a(2,:)./(1-1i*r.w); 
in.function{4}=@(o,r) real(o{4});
in.transforms={1,1,1,1};                        %%Field transforms
in.ftransforms= {1,1,1,1};                      %%Function transforms
in.olabels{1} =   '|a(\omega)|^2';
in.olabels{2} =   '|a_{in}(\omega)|^2';
in.olabels{3} =   '|a_{out}(\omega)|^2';
in.olabels{4} =   '\Re(a(\omega)-\surd{2}a_{in}(\omega)/(1+i\omega))';
in.olabels{5} =   '|a(\omega)-\surd{2}a_{in}(\omega)/(1+i\omega)|^2';
in.olabels{6} =   '|a(t)|^2';
in.compare{1}= @(w,r) 1./(1+w.^2);               %%Comparisons
in.compare{2}= @(w,r) 0.5+0*w;                   %%Comparisons
in.compare{3}= @(w,r) 0.5+0*w;                   %%Comparisons
e            = xspde(in);                        %%Stochastic toolbox
end