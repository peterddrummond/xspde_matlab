function O  =  ha(k,psi)            
%   O = ha(k,psi) Computes k-th Hadamard operator on psi
%   xSPDE functions licensed by Peter D. Drummond, (2024) - see License

szp = size(psi);
sc = [prod(szp(1:k-1)),szp(k),prod(szp(k+1:end))];
psi = reshape(psi,sc);
O   = complex(zeros(sc));
O(:,1,:) =  psi(:,1,:)+psi(:,2,:); 
O(:,2,:) =  psi(:,1,:)-psi(:,2,:); 
O   = reshape(O,szp)/sqrt(2);
end                                        %%End function call  