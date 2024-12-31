function e = SSElin2
%Uses a non-sparse SSE to solve for a linear two-mode decay

p.name       = 'SSE, N = 3,6';
p.N          = 3;
p.Om         = 1;
p.ranges     = 2;
p.nmax       = [p.N+1,2*p.N+1];
p.quantum    = 1;
p.ensembles  = [100, 10];
p.gamma{1}   = @(p)[0.5,1];
p.L{1}       = @a;
p.H          = @(psi,p) p.Om*(n(1,psi)+n(2,psi));
p.initial    = @(~,p) kron([0,0,0,1]',[0,0,0,0,0,0,1]);
p.expect{1}  = @(psi,p) n(1,psi);
p.expect{2}  = @(psi,p) n(2,psi);
p.compare{1} = @(p) p.N*exp(-p.t);
p.compare{2} = @(p) 2*p.N*exp(-2*p.t);
p.olabels    = {'n_1','n_2'};
e  = xspde(p);
end