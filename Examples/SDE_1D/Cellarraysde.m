function e = Cellarraysde()
p.fields      = {2,1}; 
p.noises      = 2;
p.initial     = {@(u,v,p) [0;2],@(u,v,p) 2;};
p.ensembles   = [100,100];
p.deriv{1}    = @(a,b,w,p) [2;1] - a + [1,0;0,0.5]*w;
p.deriv{2}    = @(a,b,w,p) - b + a(1,:);
p.observe     = {@(a,b,p) [a;b],@(a,b,p) [a.^2;b.^2]};
p.output{2}   = @(o,p)   o{2} - o{1}.^2;
p.compare{1}  = @(p) [2*(1-exp(-p.t));1+exp(-p.t);2*(1-p.t.*exp(-p.t))];
p.compare{2}  = @(p) [0.5*(1-exp(-2*p.t));0.125*(1-exp(-2*p.t));...
                      0.25*(1-(1+2*p.t+2*p.t.^2).*exp(-2*p.t))];
p.olabels{1}  = '<a(i)>,  <b>';
p.olabels{2}  = '<[\Delta a(i)]^2>,  <[\Delta b]^2>';
e = xspde(p);