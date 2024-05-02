function d = ME(rho,~,p)
%d = ME(rho,~,p) gives derivatives for master equations
%The density matrix is a sparse form of psi x psi', normalised to unity
%This is a non-stochastic version of the master equation
%The n-th dissipative operator for mode m is L_nm = sqrt(g_nm)L{n}{m}, 
%with amplitude decay rate g_nm = gamma{n}(m), and overall master equation:
%
% d_rho/dt = -i[H0,rho]+sum_nm([L_nm.rho,L_nm^dg]+[L_nm,rho.L_nm^dg])
%
%   Input:  density matrix rho, data structure p 
%   Output: derivative d
%   Needs: Hamiltonian p.H, indices p.nl, rates p.gamma, operators p.L
%   Called by: p.method
%   xSPDE functions are licensed by Peter D. Drummond, (2023) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d  = -1i*(p.H(p)*rho-rho*p.H(p));
for n = p.nl
  gm   = p.gamma{n}(p);
  szg = size(gm);
  for m = 1:numel(gm)
    g = gm(m);
    mv = m;
    if szg(1)> 1
        [mv(1),mv(2)] = ind2sub(szg,m);
    end
    if g > 0
      L    = sqrt(g)*p.L{n}(mv,p);
      LdL  = L'*L;
      d    = d + 2*(L*rho)*L' - LdL*rho - rho*LdL;
    end
  end
end
end