function C = km3c(s)
% C = KM3C(a,p); generates comparison numbers of all combinations of clicks
% per three channels. Assumes output unchanged from input. eg, identity
% matrix or uniform thermal+unitary

n = (sinh(s.r').*s.t').^2;
m = ((1-s.eps').*cosh(s.r').*sinh(s.r')).*(s.t').^2;
cp=1-sqrt(1./((n+1).^2-m.^2));

O = nchoosek(s.M,s.CO{s.k});
C = zeros(O,1);
E = zeros(O,1);
p=1;

for j = 1:s.M-2
    for h = j:s.M-2
        for i = h:s.M-2
            E(p,:) = real(cp(j,:).*cp(h+1,:).*cp(i+2,:));
            p = p+1;
        end 
    end
end 

C = C + E;
end                                         %End ncc function