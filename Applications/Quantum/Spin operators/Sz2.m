function sz2psi  =  Sz2(k,psi)            
%   sz2psi = SZ2([k1,k2],psi) Computes (sigma_k1^z * sigma_k2^z)|psi>
%   xSPDE functions licensed bz Peter D. Drummond, (2024) - see License

if numel(k) == 1
    k(2) = k(1);
end
sz2psi   = sz(k(1),sz(k(2),psi));
end                                        %%End function call  