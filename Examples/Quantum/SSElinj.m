function [e] = SSElinj
%Uses an SSE to solve for a linear decay

p.name       = 'SSE jump linear one-mode decay, N=6 initial photons';
p.nmax       = 7;
p.ranges     = 2;
p.quantum    = 1;
p.jump       = 1;
p.ensembles  = [100,1,10];
p.gamma{1}   = @(p) 0.25;
p.compare{1} = @(p) 6*exp(-0.5*p.t);
p.L{1}       = @a;
p.H          = @(psi,p) n(1,psi);
p.diffplot   = {1,1};
p.initial    = @(w,p) [0,0,0,0,0,0,1]';
p.expect{1}  = @(psi,p) n(1,psi);
p.olabels    = {'\langle N \rangle','\langle N \rangle'};
e            = xspde(p);
end