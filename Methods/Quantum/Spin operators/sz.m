function szpsi  =  sz(k,psi)            
%   szpsi = sz(k,psi) Computes k-th sigma^z operator on psi
%   xSPDE functions licensed by Peter D. Drummond, (2023) - see License

szp = size(psi); k = abs(k);
sc = [prod(szp(1:k-1)),szp(k),prod(szp(k+1:end))];
psi = reshape(psi,sc);
szpsi   = complex(zeros(sc));
szpsi(:,1,:) =  psi(:,1,:); 
szpsi(:,2,:) = -psi(:,2,:);  
szpsi   = reshape(szpsi,szp);
end                                        %%End function call  