function sypsi  =  Sy(k,psi)            
%   sypsi = SY(k,psi) Computes k-th sigma^y operator on psi
%   xSPDE functions licensed by Peter D. Drummond, (2024) - see License

szp = size(psi); k = abs(k);
sc = [prod(szp(1:k-1)),szp(k),prod(szp(k+1:end))];
psi = reshape(psi,sc);
sypsi   = complex(zeros(sc));
sypsi(:,1,:) = -1i*psi(:,2,:); 
sypsi(:,2,:) =  1i*psi(:,1,:);  
sypsi   = reshape(sypsi ,szp);
end                                        %%End function call  