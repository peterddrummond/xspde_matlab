function e  =  Peregrine()
%   e  =  Peregrine() tests xSPDE for a nonlinear Schroedinger equn.
%   Tests a partial differential equation on a finite interval for:
%   Using NN,DD,DN,ND boundary values with a spectral method
%   Tests use time dependent boundary values for a peregrine solution
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License

p.dimensions = 2;
p.noises    =  1;
p.fields    =  4;
p.order     =  2;
p.ranges    =  [10,10];
p.origins   =  [-5,-5];
p.points    =  [51,161];
p.method    =  @MP;
p.olabels   =  {'|a|^2 , DD','|a|^2 , NN','|a|^2 , DN','|a|^2 , ND'};
p.boundaries{2} = [1,1;-1,-1;1,-1;-1,1];
p.boundfun  =  @boundval;
sol         =  @(p) abs(per(p.x,p.t).^2);
p.initial   =  @(~,p) per(p.x,p.origins(1))+zeros(4,1,1);  
p.compare   =  {@(p) sol(p),@(p) sol(p),@(p) sol(p),@(p) sol(p)};
p.observe   =  {@(a,p) a(1,:),@(a,p) a(2,:),@(a,p) a(3,:),@(a,p) a(4,:)};
p.output    =  {@(o,p) abs(o{1}).^2,@(o,p) abs(o{2}).^2,...
                @(o,p) abs(o{3}).^2,@(o,p) abs(o{4}).^2};
p.name     =  'PerSpecBval';
p.steps    =  20;
p.deriv    =  @(a,w,p) 1i*a.*((conj(a).*a));
p.linear   =  @(p) 0.5*1i*p.Dx.^2;
e = xspde(p);
end                            %%end of function

function [p,varargout]  =  per(x,t)
%   Generates peregrine solutions with alpha = 1/2, beta = A0 = 1
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License

p  = exp(1i*t).*(4*(1+2*1i*t)./(1+4.*(t.^2+x.^2))-1);
if nargout  == 2
 dp = -8*x.*exp(1i*t).*(4*(1+2*1i*t)./(1+4.*(t.^2+x.^2)).^2);
 varargout{1} = dp;
end
end

function bound =  boundval(~,~,~,p)
%   Generates nonzero, time dependent boundary values
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License
[p,dp] = per(p.origins(2),p.t);
bound = {p,p;dp,-dp;p,-dp;dp,p};
end



