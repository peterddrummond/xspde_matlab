function e = BoseHubbard
%Uses a non-sparse functional method
%Includes coupling in a 1D B-H model
%Initial number state, no loss

p.name       = '4-mode linear B-H, N = [3,1,2,0]';
p.N          = [3,1,2,0];
p.modes      = 4;
p.ranges     = 5;
p.steps      = 5;
p.Om         = 1;
p.nmax       = (sum(p.N)+1)*ones(1,p.modes);
p.quantum    = 1;
p.ensembles =  [10,1,12];
p.H          = @(psi,p) BoseHubbardH(p.Om,1,psi,p);
p.initial    = @(~,p) Mknumber(p.N,p);
p.expect{1}  = @(psi,p) N(1,psi);
p.expect{2}  = @(psi,p) N(2,psi);
p.expect{3}  = @(psi,p) N(3,psi);
p.expect{4}  = @(psi,p) N(4,psi);
p.olabels    = {'< n_1 >','< n_2 >','< n_3 >','< n_4 >'};
e  = xspde(p);
end