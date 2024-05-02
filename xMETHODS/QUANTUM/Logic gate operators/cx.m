function O  =  cx(k,psi)            
%   O = cx(k,psi) Computes (k1,k2)-th controlled not operator on psi
%   xSPDE functions licensed by Peter D. Drummond, (2024) - see License

szp = size(psi);
if numel(k) == 1 || k(1)==k(2)
    error('CNOT gate has identical indices %d\n',k(1));
end
j = k(2);
k = k(1);
if k<j
  sc = [prod(szp(1:k-1)),2,prod(szp(k+1:j-1)),2,prod(szp(j+1:end))];
  psi = reshape(psi,sc);
  O   = complex(zeros(sc));
  O(:,2,:,1,:) =  psi(:,2,:,2,:); 
  O(:,2,:,2,:) =  psi(:,2,:,1,:); 
else
  sc = [prod(szp(1:j-1)),2,prod(szp(j+1:k-1)),2,prod(szp(k+1:end))];
  psi = reshape(psi,sc);
  O   = complex(zeros(sc));
  O(:,1,:,2,:) =  psi(:,2,:,2,:); 
  O(:,2,:,2,:) =  psi(:,1,:,2,:);
end
O   = reshape(O,szp);
end                                        %%End function call  