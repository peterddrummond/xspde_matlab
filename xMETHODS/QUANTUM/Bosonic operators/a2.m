function a2psi  =  a2(k,psi)            
%   a2psi = a2(k,psi) is a quantum operator function
%   for bosonic annihilation/creation operators squared 
%   Returns  a_k^2|psi> if k>0 
%   Returns  a_|k|^{dag 2}|psi> if k<0 (hermitian conjugate)
%   Licensed by Peter D. Drummond, (2024) - see License

szp  = size(psi);
km   = abs(k); 
sz1  = [prod(szp(1:km-1)),szp(km),prod(szp(km+1:end))];
a2psi = reshape(psi,sz1);
if k  > 0
  for j = 1:(szp(k)-2)
    a2psi(:,j,:) = sqrt(j*(j+1))*a2psi(:,j+2,:);   
  end
  a2psi(:,szp(k),:)   = 0;
  a2psi(:,szp(k)-1,:) = 0;
else
  for j = 1:(szp(km)-2)
     a2psi(:,j+2,:) = sqrt(j*(j+1))*a2psi(:,j,:);
   end
   a2psi(:,1,:) = 0;
   a2psi(:,2,:) = 0;
end  
a2psi = reshape(a2psi,szp);
end                                              %%End function call  