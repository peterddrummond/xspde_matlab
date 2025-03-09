function C = Km(a,p)
% C = KM(a,p) calculates a click probability product moment
% This is plotted with increasing correlation order, from O(1),..O(i)..
% Hence, C(i) = <k(1).k(2)...k(O(i))>
% Here O is an input vector of orders plotted

% INITIALIZE PROJECTORS

if p.phase ~= 1                                 %check method is OK     
    error('qSIM supports click probabilities only for +P method = 1');
end                                              %end check method is OK
O = p.correl{p.nobserve};                             %store moment vector    
Ol = length(O);                                  %store output length                 
n = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);    %get normal ordered number
k  = 1-exp(-n);                                  %click projector

% COMPUTE CORRELATIONS

C = zeros(Ol,p.ensembles(1));                    %loop over output index
for i = 1:Ol                                     %loop over output index
    C(i,:) = real(prod(k(1:O(i),:),1));          %get mmoment of k
end                                              %end loop over index
end                                              %end function