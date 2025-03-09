function C = Kmsub(a,p)

% C = KMSUB(a,p); Generates all a subset of all possible combinations of
% output mode click correlations for moments of sequential modes. This is
% done for graphical simplificity. 
% Example, for a second order correlation, the first data point is a 
% product of clicks from modes 1, 2, the second data point is a product 
% of clicks from modes 2,3, third data point is a product of clicks from 
% modes 3,4 etc.

% INITIALIZE PROJECTORS

if p.phase ~= 1                                 %check method is OK     
    error('qSIM supports click probabilities only for +P method = 1');
end                                              %end check method is OK
               
n = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);                %get normal ordered number
k  = 1-exp(-n);

% Compute correlations

Ord = p.CO{p.nobserve} - 1;
NM = p.modes - Ord;

C = zeros(NM, p.ensembles(1));

for i = 1:NM
    C(i,:) = real(prod(k(i:i+Ord,:),1));
end 
end