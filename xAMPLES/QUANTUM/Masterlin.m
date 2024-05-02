function [e] = masterlin
%Uses a master equation to solve for a linear decay

p.name       = 'master equation linear one-mode decay, N=6 initial photons';
p.N          = 6;
p.nmax       = p.N+1;
p.modes      = 1;
p.ranges     = 2;
p.quantum    = 2;
p.gamma{1}   = @(p) 0.25;
p.compare{1} = @(p) p.N*exp(-2*p.t*p.gamma{1});
p.L{1}       = p.a;
p.H          = @(p) p.a{1}'*p.a{1};
p.diffplot   = {1,1};
p.initial    = @(w,p) kron([0;0;0;0;0;0;1],[0,0,0,0,0,0,1]);
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.olabels    = {'\langle n \rangle'};
e            = xspde(p);
end