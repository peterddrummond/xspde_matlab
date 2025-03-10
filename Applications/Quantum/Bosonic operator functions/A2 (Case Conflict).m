function a2psi  =  A2(k,psi)            
%   a2psi = A2(k,psi) is a quantum operator function
%   for bosonic annihilation/creation operators squared 
%   Returns  a_k^2|psi> if k>0 
%   Returns  a_|k|^{dag 2}|psi> if k<0 (hermitian conjugate)
%   Licensed by Peter D. Drummond, (2024) - see License

a2psi = A(k,A(k,psi));
end                                              %%End function call  