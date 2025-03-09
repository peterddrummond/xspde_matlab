function [e] = SSElin
%Uses an SSE to solve for a linear decay

p.name       = 'SSE linear decay, N=6';
p.nmax       = 7;
p.ranges     = 2;
p.quantum    = 1;
p.ensembles  = [100,1,10];
p.gamma{1}   = @(p) 0.25;
p.compare{1} = @(p) 6*exp(-2*p.t*p.gamma{1}(p));
p.L{1}       = @A;
p.H          = @(psi,p) N(1,psi);
p.diffplot   = {1,1};
p.initial    = @(w,p) Mknumber(6,p);
p.expect{1}  = @(psi,p) N(1,psi);
p.olabels    = {'\langle N \rangle','\langle N \rangle'};
e            = xspde(p);
end