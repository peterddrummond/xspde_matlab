function e = Masterlin2
%Uses a master equation to solve for a linear two-mode decay

p.name       = 'Master equation, N = 3,6';
p.N          = [3,6];
p.noises     = 4;
p.ranges     = 2;
p.nmax       = p.N+1;
p.a          = Mkbose(p);
p.quantum    = 2;
p.gamma      = @(p)[0.5,1];
p.L{1}       = @(m,p) p.a{m};
p.initial    = @(~,p) Mknumber(p.N,p);
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.expect{2}  = @(p) p.a{2}'*p.a{2};
p.compare{1} = @(p) p.N(1)*exp(-p.t);
p.compare{2} = @(p) p.N(2)*exp(-2*p.t);
p.diffplot   = {1,1};
p.olabels    = {'< n_1 >','< n_2 >'};
e  = xspde(p);
end