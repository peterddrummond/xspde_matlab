function [e] = SSElin
%Uses an SSE to solve for a linear decay

p.name       = 'SSE linear decay, N=6';
p.nmax       = 7;
p.ranges     = 2;
p.quantum    = 1;
p.ensembles  = [100,1,10];
p.gamma{1}   = @(p) 0.25;
p.compare{1} = @(p) 6*exp(-2*p.t*p.gamma{1}(p));
p.L{1}       = @a;
p.H          = @(psi,p) n(1,psi);
p.diffplot   = {1,1};
p.initial    = @(w,p) [0,0,0,0,0,0,1]';
p.expect{1}  = @(psi,p) n(1,psi);
p.olabels    = {'\langle N \rangle','\langle N \rangle'};
e            = xspde(p);
end