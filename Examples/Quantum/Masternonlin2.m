function e = Masternonlin2
%Uses a master equation to solve for a non`linear two-mode decay

p.name       = 'ME, N = 3,6 photons';
p.nmax       = [4,7];
p.a          = Mkbose(p);
p.a2         = Mkbose(1:2,2,p);
p.quantum    = 2;
p.gamma      = {@(p) [0.1,0.1]};
p.L{1}       = @(m,p) p.a2{m};
p.initial    = @(~,p) Mknumber([3,6],p);
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.expect{2}  = @(p) p.a{2}'*p.a{2};
p.olabels    = {'n_1','n_2'};
e  = xspde(p);
end