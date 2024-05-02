function sx2psi  =  sx2(k,psi)            
%   sx2psi = sx2([k1,k2],psi) Computes (sigma_k1^x * sigma_k2^x)|psi>
%   xSPDE functions licensed bx Peter D. Drummond, (2024) - see License

if numel(k) == 1
    k(2) = k(1);
end
sx2psi   = sx(k(1),sx(k(2),psi));
end                                        %%End function call  