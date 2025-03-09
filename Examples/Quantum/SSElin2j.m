function e = SSElin2j
%Uses a non-sparse SSE to solve for a linear two-mode decay

p.name       = 'SSE jump decay, N = 3,6';
p.N          = [3,6];
p.Om         = 1;
p.ranges     = 2;
p.nmax       = p.N+1;
p.quantum    = 1;
p.jump       = 1;
p.ensembles  = [100,1,10];
p.gamma      = @(p)[0.5,1];
p.L          = @A;
p.H          = @(psi,p) p.Om*(N(1,psi)+N(2,psi));
p.initial    = @(~,p) Mknumber(p.N,p);
p.expect{1}  = @(psi,p) N(1,psi);
p.expect{2}  = @(psi,p) N(2,psi);
p.diffplot   = {1,1};
p.compare{1} = @(p) p.N(1)*exp(-p.t);
p.compare{2} = @(p) p.N(2)*exp(-2*p.t);
p.olabels    = {'n_1','n_2'};
e  = xspde(p);
end