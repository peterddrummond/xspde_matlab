function C = n1lsc(p)
% C = N1LSC(p); generates comparison count probabilities in one partition
% of length p.modes, with squeezed inputs,  p.thermal=0, and unitary transmission
if p.thermal ~= 0
    error('n1lsc uses uniform squeezed state inputs');
end


M = p.modes*0.5;                    %Modes/2
t = p.tr^2;                         %Squared amplitude loss
p1  = 1/(1+mean((sinh(p.sqz)).^2,2)); %probability one 
lp2 = log(1-p1);                    %probability two

C = zeros(p.max+1,1);               %Initialize total count vector
OC = C;                             %Initialize odd count vector
opC = zeros(p.max+1,1);
opC(1,1) = log(M) + lp2;
EC = C;
epC = M*log(p1)+zeros(p.max+1,1);

%loop over even counts. Here, we use \tilde{m} = m/2
for m = 1:p.max
    ehg = log(hypergeom([m + 0.5, M + m],0.5,(1-t)^2.*(1-p1)));
    t2 = 2*m*log(t);
    epC(m+1,1) = epC(m,1)+log(1+(M-1)/m)+lp2;
    EC(m+1,1) = t2 + ehg + epC(m+1,1);
end
EC(1,1) = epC(1,1) + log(hypergeom([0.5,M],0.5,(1-t)^2.*(1-p1)));%Replace add hypergeometric to m=0 count
EC = exp(EC);


%loop over odd counts. Here, we use m' = (m+1)/2
for m = 1:p.max+1
    % ohg = log(hypergeom([0.5*(m + 2), M + 0.5*(m + 1)],1.5,(1-t)^2.*(1-p1)));
    % tm = m*log(t);
    ohg = log(hypergeom([m + 0.5, M + m],1.5,(1-t)^2.*(1-p1)));
    tm = (2*m - 1)*log(t);
    if m ==1 
        opC(m,1) = log(M) + lp2;
    else 
        opC(m,1) = opC(m-1,1) + log(1+(M-1)/m)+lp2;
    end 
    OC(m,1) = M*log(p1) + log(2*m) + tm + ohg + opC(m,1);
end 
OC = (1-t).*exp(OC);

j = 1;
k = 1;
for i = 1:p.max+1
    if mod(i,2) == 1
        C(i,1) = EC(j,1);
        j = j+1;
    elseif mod(i,2) == 0
        C(i,1) = OC(k,1);
        k = k+1;
    end 
end 

