function npsi  =  n(k,psi)            
%   npsi = n(k,psi) is a quantum operator function
%   for bosonic number operators
%   Returns  n_k|psi> = a_k^dag*a_k|psi> for scalar k
%   Returns  a_i^dag*a_j|psi> for vector k =  [i,j] 
%   Returns  a_j^dag*a_i|psi> for vector k = -[i,j] 
%   Negative first mode number ==> conjugate operator
%   Licensed by Peter D. Drummond, (2024) - see License

nk = numel(k);
if nk == 1
  k = abs(k(1));
  sz = size(psi);
  sc = [prod(sz(1:k-1)),sz(k),prod(sz(k+1:end))];
  psi = reshape(psi,sc);
  npsi   = zeros(sc);
  for n = 0:(sz(k)-1)
    npsi(:,n+1,:) = n*psi(:,n+1,:);
  end
  npsi = reshape(npsi,sz);
else
  if k(1) > 0
    npsi = a(-k(1),a(k(2),psi));
  else
    npsi = a(k(2),a(-k(1),psi));
  end
end
end                                        %%End function call  