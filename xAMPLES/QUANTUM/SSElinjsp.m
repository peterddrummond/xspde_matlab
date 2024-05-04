function [e] = SSElinjsp
%Uses an SSE to solve for a linear decay with sparse methods

p.name       = 'SSE sparse jump decay, N=6';
p.N          = 6;
p.nmax       = p.N+1;
p.ranges     = 2;
p.a          = mkbose(p);
p.quantum    = 1;
p.sparse     = 1;
p.jump       = 1;
p.ensembles  = [100,1,10];
p.gamma{1}   = @(p) 0.25;
p.compare{1} = @(p) p.N*exp(-2*p.t*p.gamma{1}(p));
p.L{1}       = @(m,p) p.a{m};
p.H          = @(p) p.a{1}'*p.a{1};
p.diffplot   = {1,1};
p.initial    = @(w,p) (p.a{1}')^6*mkvacs(p)/sqrt(720);
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.olabels    = {'\langle N \rangle','\langle N \rangle'};
e            = xspde(p);
end