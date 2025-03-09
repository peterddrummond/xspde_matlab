function sy2psi  =  Sy2(k,psi)            
%   sy2psi = SY2([k1,k2],psi) Computes (sigma_k1^y * sigma_k2^y)|psi>
%   xSPDE functions licensed by Peter D. Drummond, (2024) - see License

if numel(k) == 1
    k(2) = k(1);
end
sy2psi   = sy(k(1),sy(k(2),psi));
end                                        %%End function call  