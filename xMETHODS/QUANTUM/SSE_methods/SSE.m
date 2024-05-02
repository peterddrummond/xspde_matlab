function d = SSE(psi,w,p)
%d = SSE(psi,w,p) gives derivatives for a stochastic Schroedinger equation
%The wavefunction is stored as psi(n_1,n_2,n_k...), normalised to unity
%The n-th dissipative operator for mode k is p.L{n}(psi,[k1,..kx],p)
%The n-th amplitude decay rate is p.gamma{n}(k), where k is a mode index.
%If required, k=[k1,..kx] is a vector index for more than one mode.
%This index is passed to the reservoir coupling operator p.L{n}(psi,k)
%The non-dissipative Hamiltonian acting on psi is p.H(psi,p)
%The equation is in the Stratonovich form of stochastic calculus.
%
%Master equn: d rho/dt=-i[H,rho]+sum_n([L_n.rho,L_n^dg]+[L_n,rho.L_n^dg])
%Where: L_n = sum_k(gamma_n(k).L_n(k)
%
%SSE: d psi/dt = -iH. psi + sum_n [(L_n - <L_n>).(z + 2<L_n^dg>)
%                                 + <L_n^dg L_n>-L_n^dg L_n].psi
%where <<dz dz^*>> = 2dt
%
%See: Diosi et al, PRA58, 1699, (1998) and others, note scaling of gamma.
%For best results, combine this with a projected algorithm, p.MPS
%   Input:  wavefunction psi, real noises w, data structure p 
%   Output: derivative d
%   Called by: p.method
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l = 1;                                           %%Initialize noise index
nd = p.nfields;                                  %%Store number of modes
if p.sparse                                      %%If using sparse methods
  d  = -1i*p.H(p)*psi;                           %%Multiply psi by H
else                                             %%Else not sparse methods
  d  = -1i*p.H(psi,p);                           %%Use H.psi as a function
end                                              %%End if sparse methods
for n = p.nl                                     %%Loop over loss terms
  gm   = p.gamma{n}(p);                          %%Get n-th loss rate
  szg = size(gm);                                %%Get size of loss rate
  for m = 1:numel(gm)                            %%Loop on n-th loss vector
    g = gm(m);                                   %%Store m-th loss rate
    theta = p.theta{n}(m);                       %%Store m-th loss phase
    mv = m;
    if szg(1)> 1
        [mv(1),mv(2)] = ind2sub(szg,m);
    end
    if g > 0   
      if theta == 0
        z  = sqrt(g)*(w(l,:)+1i*w(l+1,:));
        l  = l + 2;
      else
        z  = sqrt(2*g)*theta*w(l,:);
        l  = l + 1;
      end
      if p.sparse
        L  = p.L{n}(mv,p);
        X  = L' + theta^2*L;
        L  = L*psi;
        EX = sum(conj(psi).*(X*psi));
        XL = X*L;
      else
        z  = reshape(z,[ones(1,nd),p.ensembles(1)]);
        L  = p.L{n}(mv,psi);
        X  = p.L{n}(-mv,psi) + theta^2*L;
        XL = p.L{n}(-mv,L) + theta^2*p.L{n}(mv,L);
        EX = sum(conj(psi).*X,1:nd);
      end
      d  = d - g*(XL-sum(conj(psi).*XL,1:nd).*psi);
      EL = sum(conj(psi).*L,1:nd);
      L  = L - EL.*psi;   
      d  = d + (2*g*EX+z).*L ;     
    end
  end
end
end