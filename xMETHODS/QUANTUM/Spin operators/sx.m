function sxpsi =  sx(k,psi)            
%   sxpsi = sx(k,psi) Computes k-th sigma^x operator on psi
%   xSPDE functions licensed by Peter D. Drummond, (2024) - see License

szp = size(psi); k = abs(k);
sc = [prod(szp(1:k-1)),szp(k),prod(szp(k+1:end))];
psi = reshape(psi,sc);
sxpsi   = complex(zeros(sc));
sxpsi(:,1,:) = psi(:,2,:); 
sxpsi(:,2,:) = psi(:,1,:);  
sxpsi   = reshape(sxpsi,szp);
end                                        %%End function call  