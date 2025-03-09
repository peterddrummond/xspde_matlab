function O  =  P8(k,psi)            
%   O = P8(k,psi) Computes k-th pi/8 operator on psi
%   xSPDE functions licensed by Peter D. Drummond, (2024) - see License

szp = size(psi);
sc = [prod(szp(1:k-1)),szp(k),prod(szp(k+1:end))];
psi = reshape(psi,sc);
O   = complex(zeros(sc));
O(:,1,:) =  psi(:,1,:); 
O(:,2,:) = exp(1i*pi/4)*psi(:,2,:);  
O   = reshape(O,szp);
end                                        %%End function call  