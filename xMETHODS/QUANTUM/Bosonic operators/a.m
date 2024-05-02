function apsi  =  a(k,psi)            
%   apsi = a(k,psi) is a quantum operator function
%   for bosonic annihilation/creation operators
%   Returns  a_k|psi> if k>0 
%   Returns  a_|k|^dag|psi> if k<0 (hermitian conjugate)
%   Licensed by Peter D. Drummond, (2024) - see License

szp  = size(psi);
km   = abs(k);
sz1  = [prod(szp(1:km-1)),szp(km),prod(szp(km+1:end))];
apsi = reshape(psi,sz1);
n    = reshape((0:(szp(km)-1)),[1,szp(km),1]);
if k > 0
    apsi = circshift(sqrt(n).*apsi,[0,-1,0]);
else
    apsi = sqrt(n).*circshift(apsi,[0,1,0]);
end
apsi = reshape(apsi,szp);
end                                              %%End function call  