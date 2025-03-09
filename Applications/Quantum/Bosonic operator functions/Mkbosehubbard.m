function h = MkBoseHubbard(O,K,p)
% H = MkBOSEHUBBARD(p) gives the Bose-Hubbard  Hamiltonian
% in a sparse form. 
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
  for m = 1:p.modes
    h = h+K*(p.a{m}'*p.a{m})^2;
  end
else
  for m = 1:p.modes
    h = h+K(m)*(p.a{m}'*p.a{m})^2;
  end
end
if isscalar(O)
  for m = 1:p.modes-1
      h = h+O*(p.a{m}'*p.a{m+1}+p.a{m+1}'*p.a{m})/2;
  end
elseif isvector(O)
  for m = 1:max(p.modes,length(O))
    if m < p.modes
      h = h+O*(p.a{m}'*p.a{m+1}+p.a{m+1}'*p.a{m})/2;
    else
      h = h+O*(p.a{m}'*p.a{1}+p.a{1}'*p.a{m})/2;
    end
  end
else
  for m = 1:p.modes
    for n = 1:p.modes
      h = h+O(m,n)*p.a{m}'*p.a{n}/2;
    end
  end
end
end