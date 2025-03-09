function h = BoseHubbardH(O,K,psi,p)
% H = BOSEHUBBARDH(O,K,psi,p) gives the Bose-Hubbard Hamiltonian
% in a functional form. 
%
% H = (1/2)\sum _nm O_mn a'(m)a(n) + \sum _n K_n [a'(n)a(n)]^2 
%
% Scalar nonlinear parameters K are extended to vectors
% Scalar couplings O are extended to vectors of length modes-1
% Vector couplings O are assumed to hold for linear nearest neighbors
% Output: Hamiltonian H
% Needs: Couplings O, K
%
% xSPDE functions are licensed by Peter D. Drummond, (2025) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=0;
if isscalar(K)
    K = K*ones(1,p.modes);
end
for m = 1:p.modes
    h = h+K(m)*N(m,N(m,psi));
end
if isscalar(O)
    O = O*ones(1,p.modes-1);
end
if isvector(O)
  for m = 1:min(p.modes,length(O))
    if m < p.modes
     h = h+O(m)*(A(-m,A(1+m,psi))+A(-(m+1),A(m,psi)))/2;
    else
      h = h+O(m)*(A(-m,A(1,psi))+A(-1,A(m,psi)))/2;
    end
  end
else
  for m = 1:p.modes
    for n1 = 1:p.modes
      h = h+O(m,n1)*A(-m,A(n1,psi))/2;
    end
  end
end
end