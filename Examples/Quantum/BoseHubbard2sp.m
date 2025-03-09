function e = BoseHubbard2sp
%Uses a Schroedinger equation with an initial number state
%example from Carusotto,Castin,Dalibard,arXiv/0003399
%Includes coupling in a B-H model
p.name       = '2-mode coupled B-H, N = 17';
p.N          = [17,0];
p.modes      = 2;
p.points     = 501;
p.order      = 4;
p.K          = 0.1;
p.O          = 1.0;
p.nmax       = (sum(p.N)+1)*ones(1,2);
p.a          = Mkbose(p);
p.BH         = Mkbosehubbard(p.O,p.K,p);
p.quantum    = 1;
p.sparse     = 1;
p.H          = @(p) p.BH;
p.initial    = @(~,p) Mknumber(p.N,p);
p.expect{1}  = @(p) p.a{1}'*p.a{1}/p.N(1);
p.expect{2}  = @(p) p.a{2}'*p.a{2}/p.N(1);
p.olabels    = {'< n_1 >/N','< n_2 >/N'};
e  = xspde(p);

end