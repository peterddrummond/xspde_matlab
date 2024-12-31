function C = km3(a,p)
% C = km3(a,p); generates ALL third-order click correlation moments for
% sequential modes where j<k<h. 


% INITIALIZE PROJECTORS

if p.phase ~= 1
    error('Qsim supports click probabilities only for +P method = 1');
end

n = a(1:p.modes,:).*a(p.modes+1:2*p.modes,:);                %normal ordered numbers
k  = 1-exp(-n);                                  %click projector

% COMPUTE CORRELATIONS

O = nchoosek(p.modes,p.CO{p.nobserve});
C = zeros(O,p.ensembles(1));
E = zeros(O,p.ensembles(1));
s=1;

for j = 1:p.modes-2
    for h = j:p.modes-2
        for i = h:p.modes-2
            E(s,:) = real(k(j,:).*k(h+1,:).*k(i+2,:));
            s = s+1;
        end 
    end
end 

C = C + E;
end 