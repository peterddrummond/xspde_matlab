function [e] = Masterlin
%Uses a master equation to solve for a linear decay

p.name       = 'Master equation, N=6';
p.N          = 6;
p.nmax       = p.N+1;
p.modes      = 1;
p.ranges     = 2;
p.quantum    = 2;
p.a          = mkbose(p);
p.gamma      = @(p) 0.25;
p.compare{1} = @(p) p.N*exp(-p.t*0.5);
p.L{1}       = @(m,p) p.a{m};
p.H          = @(p) p.a{1}'*p.a{1};
p.diffplot   = {1,1};
p.initial    = @(w,p) kron([0;0;0;0;0;0;1],[0,0,0,0,0,0,1]);
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.olabels    = {'\langle n \rangle'};
p.rawdata = 1;
e            = xspde(p);
end