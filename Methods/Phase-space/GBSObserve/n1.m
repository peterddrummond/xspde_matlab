function C = n1(a,p)
% C = n1(a,p); generates count probabilities in a single partition
% return probability

if p.phase ~= 1
    error('n1 requires +P method = 1');
end
np = sum(a(1:p.modes,:).*a(p.modes+1:2*p.modes,:),1);        %photon numbers
%C = zeros(p.max+1,1);                            %total count
lp = log(np);
for m = 1:p.max+1
    %C(m,1) = mean(real(exp(-np)),2);
    C(m,:) = real(exp(-np));
    np = np+log(m)-lp;
end                                              %End n1 function