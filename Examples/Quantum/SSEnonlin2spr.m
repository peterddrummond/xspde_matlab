function e = SSEnonlin2spr
%Uses a master equation to solve for a nonlinear two-mode decay

p.name       = 'SSE real sparse nonlin, N = 3,6';
p.nmax       = [4,7];
p.steps      = 4;
p.a          = mkbose(p);
p.a2         = mkbose(1:2,2,p);
p.ensembles  = [10,10,10];
p.quantum    = 1;
p.sparse     = 1;
p.theta     = [1,1];
p.gamma      = {@(p) [0.01,0.01],@(p)[.5,.25]};
p.L          = {@(m,p) p.a{m},@(m,p) p.a2{m}};
p.initial    = @(~,p) kron([0,0,0,1],[0,0,0,0,0,0,1])';
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.expect{2}  = @(p) p.a{2}'*p.a{2};
p.olabels    = {'n_1','n_2'};
e  = xspde(p);
end