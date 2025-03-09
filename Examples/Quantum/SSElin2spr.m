function e = SSElin2spr
%Uses a sparse SSE to solve for a linear two-mode decay

p.name       = 'SSE sparse real, N = 3,6';
p.N          = [3,6];
p.Om         = 1;
p.noises     = 4;
p.ranges     = 2;
p.nmax       = p.N+1;
p.a          = Mkbose(p);
p.Hmat       = p.Om*(p.a{1}'*p.a{1}+p.a{2}'*p.a{2});
p.quantum    = 1;
p.sparse     = 1;
p.ensembles  = [100,1,10];
p.theta{1}   = [1,1];
p.gamma{1}   = @(p) [0.5,1]*p.t;
p.L{1}       = @(m,p) p.a{m};
p.H          = @(p) p.Hmat;
p.initial    = @(~,p) Mknumber(p.N,p);
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.expect{2}  = @(p) p.a{2}'*p.a{2};
p.compare{1} = @(p) p.N(1)*exp(-p.t.^2/2);
p.compare{2} = @(p) p.N(2)*exp(-p.t.^2);
p.diffplot   = {1,1};
p.olabels    = {'< n_1 >','< n_2 >'};
e  = xspde(p);
end