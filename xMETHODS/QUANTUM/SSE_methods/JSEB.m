function a = SSEJumpB(a,xi,p)
%a = SSEJUMPB(a,xi,p) gives a jump in a stochastic Schroedinger equation
% where this is the stochastic jump term
%The n-th dissipative operator for mode m is L_nm = sqrt(g_nm)L{n}{m}, 
%with amplitude decay rate g_nm = gamma{n}(m), and overall master equation:
%
% d_rho/dt = -i[H0,rho]+sum_nm([L_nm.rho,L_nm^dg]+[L_nm,rho.L_nm^dg])
%
% The equivalent stochastic Schroedinger equation or SSE is:
%
% d_psi = (JumpA.dt)+sum_nm [L_nm/<L_nm^dg L_nm> -1]dN_nm 
% 
% with:  <dN_nm^2> = 2dt<L_nm^dg L_nm>, dN_nm = 0 or 1
%  
%  See: Gardiner and Zoller and others, note gamma is amplitude decay.
%  If required, the last psi index is an ensemble index
%  Should be used with a projected algorithm, MPS or RK4S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Input:  wavefunction psi, uniform noise w=[0,1/dt], data structure p 
%   Output: new wavefunction psi
%   Called by: p.method
%   xSPDE functions are licensed by Peter D. Drummond, (2023) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = a{1};
d = zeros(size(psi));
eps = 1.e-20;                                    %%prevents divide by zero
l = 0;
sz = [ones(1,p.nfields),p.ensembles(1)];
for n = p.nl
  gm   = p.gamma{n}(p);
  szg = size(gm);
  for m = 1:numel(gm)
    g = gm(m);
    mv = m;
    if szg(1)> 1
        [mv(1),mv(2)] = ind2sub(szg,m);
    end
    if g>0
      l = l+1;
      if p.sparse
        L = sqrt(g)*p.L{n}(mv,p)*psi;
      else
        L = sqrt(g)*p.L{n}(mv,psi);
      end
      ELdL = sum(conj(L).*L,1:p.nfields);
      dN   = 2*ELdL > reshape(xi{1}(l,:),sz);
      L    = L./sqrt(ELdL+eps);
      d    = d + dN.*(L-psi);
    end
  end
end
a{1} = psi + d;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%