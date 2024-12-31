function [e] = SolitonDS()
p.dimensions = 2;
p.fields     = 4;
p.method     = @MP;
p.checks     = 0;
p.points     = [101,21];
p.steps      = 20;
p.ranges     = [1,2];
x0           = sech(1); 
x1           = tanh(1)*sech(1);
p.initial    = @(v,p) sech(p.x);
p.observe    = {@(a,p) abs(a(1,:)).^2, @(a,p) abs(a(2,:)).^2,...
                @(a,p) abs(a(3,:)).^2, @(a,p) abs(a(4,:)).^2};
p.olabels    = {'|a|^2, DD','|a|^2, NN','|a|^2, DN','|a|^2, ND'};
p.boundaries = {0,[1,1;-1,-1;1,-1;-1,1]};
p.boundval   = {0,{x0,x0;x1,-x1;x0,-x1;x1,x0}} ;
p.name       = 'SolitonDS';
p.compare    = {@(p) sech(p.x).^2,@(p) sech(p.x).^2,...
                @(p) sech(p.x).^2,@(p) sech(p.x).^2};
p.deriv      = @(a,~,p) 1i*(a.*(conj(a).*a)+0.5*(DS(a,2,p)-a));
e            = xspde(p);
end