function d = JSEA(psi,~,p)
%d = JSEA(psi,w,p) gives derivatives for a stochastic Schroedinger equation
% in a jump form, where this is the non-stochastic or "A" derivative
% The n-th dissipative operator for mode m is L_nm = sqrt(g_nm)L{n}{m}, 
% amplitude decay rate g_nm = gamma{n}(m), and the master equation is:
%
% d_rho/dt = -i[H0,rho]+sum_nm([L_nm.rho,L_nm^dg]+[L_nm,rho.L_nm^dg])
%
% The equivalent stochastic Schroedinger equation or SSE is:
%
% d_psi/dt = -iH. psi + sum_nm [ <L_nm^dg L_nm>].psi +(projection)
% 
% with:  H = H0 - i sum_nm L_nm^dg L_nm
%  
%  See: Gardiner and Zoller and others, note gamma is amplitude decay.
%  If required, the last psi index is an ensemble index
%  Combine this with a projected algorithm
%   Input:  wavefunction psi,  data structure p 
%   Output: derivative d
%   Called by: p.method
%   xSPDE functions are licensed by Peter D. Drummond, (2023) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if p.sparse
  d  = -1i*p.H(p)*psi;
else
  d  = -1i*p.H(psi,p);
end
for n = p.nl
  gm  = p.gamma{n}(p);
  szg = size(gm);
  for m = 1:numel(gm)
    if gm(m) > 0
    g = gm(m);
    mv = m;
    if szg(1)> 1
        [mv(1),mv(2)] = ind2sub(szg,m);
    end
      if p.sparse
        L   = p.L{n}(mv,p);
        LdL = (L'*L)*psi;
        L   = L*psi;
      else
        L   = p.L{n}(mv,psi);
        LdL = p.L{n}(-mv,L);
      end
      d  = d + g*(sum(conj(L).*L,1:p.nfields).*psi - LdL);
    end
  end
end 
end