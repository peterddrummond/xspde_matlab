function C = km2c(s)
% C = KM2C(a,p); generates comparison numbers of all combinations of clicks
% per two channels. Assumes output unchanged from input. eg, identity
% matrix or uniform thermal+unitary

n = (sinh(s.r').*s.t').^2;
m = ((1-s.eps').*cosh(s.r').*sinh(s.r')).*(s.t').^2;
cp=1-sqrt(1./((n+1).^2-m.^2));


O = nchoosek(s.M,s.CO{s.k});
C = zeros(O,1);
A = cell(1,s.M-1);


for i=1:s.M-1
    for j=i:s.M-1
        A{i}(j-i+1,:) = real(cp(i,:).*cp(j+1,:));
    end 
end 
C = C + cat(1,A{:});
end                                         %End ncc function