function C = km2(a,p)
% C = km2(a,p); generates ALL second-order click correlation moments for
% sequential modes where j<k. 


% INITIALIZE PROJECTORS

if p.phase ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end

n = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);                %normal ordered numbers
k  = 1-exp(-n);                                  %click projector

% COMPUTE CORRELATIONS

O = nchoosek(p.modes,p.CO{p.nobserve});
C = zeros(O,p.ensembles(1));
A = cell(1,p.modes-1);


for i=1:p.modes-1
    for j=i:p.modes-1
        A{i}(j-i+1,:) = real(k(i,:).*k(j+1,:));
    end 
end 
C = C + cat(1,A{:});
end 