function e = SSEnonlin2sp
%Uses a real sparse SSE to solve for nonlinear two-mode decay 
p.name       = 'Real sparse SSE, M=2, N=3,6';
p.nmax       = [4,7];
p.steps      = 8;
p.a          = mkbose(p);
p.a2         = mkbose(1:2,2,p);
p.ensembles  = [10,10,10];
p.quantum    = 1;
p.sparse     = 1;
p.gamma      = {@(~) [0.01,0.01],@(~) [.5,.25]};
p.alpha      = {[1,1],[1,1]};
p.L          = {p.a,p.a2};
p.initial    = @(~,p) kron([0,0,0,1],[0,0,0,0,0,0,1])';
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.expect{2}  = @(p) p.a{2}'*p.a{2};
p.olabels    = {'n_1','n_2'};
e  = xspde(p);
end