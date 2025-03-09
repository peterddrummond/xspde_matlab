function e = BoseHubbardME
%Uses a master equation with an initial number state
%Includes coupling and linear and nonlinear decay in a B-H model

p.name       = '3-mode coupled B-H, with linear and nonlinear damping';
p.N          = [3,0,0];
p.modes      = 3;
p.steps      = 5;
p.O          = 1;
p.nmax       = (sum(p.N)+1)*ones(1,p.modes);
p.a          = Mkbose(p);
p.a2         = Mkbose(1:p.modes,2,p);
p.BH         = Mkbosehubbard(p.O,1,p);
p.quantum    = 2;
p.gamma      = {@(p) [0.01,0.01,0.01],@(p) [0.2,0.2,0.2]};
p.L          = {@(m,p) p.a{m},@(m,p) p.a2{m}};
p.H          = @(p) p.BH;
p.initial    = @(~,p) Mknumber(p.N,p);
p.expect{1}  = @(p) p.a{1}'*p.a{1};
p.expect{2}  = @(p) p.a{2}'*p.a{2};
p.expect{3}  = @(p) p.a{3}'*p.a{3};
p.expect{4}  = @(p) p.a{1}'*p.a{1}+p.a{2}'*p.a{2}+p.a{3}'*p.a{3};
p.olabels    = {'< n_1 >','< n_2 >','< n_3 >','< N >'};
e  = xspde(p);
end