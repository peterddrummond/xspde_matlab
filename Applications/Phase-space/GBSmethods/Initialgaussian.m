function a = Initialgaussian(varargin)
%   a = INITIALGAUSSIAN(w,p)  generates phase-space initial gaussians
%   Input:  'w' is the noise, if needed,
%           'p' is the structure of parameters used
%   Output: 'a' is the phase-space output, [alp;(bet);(W)] OR [n;(W)]
% First index: up to p.modes or 2*p.modes
% Second  index goes up to p.ensembles(1).  
% First the input amplitude p.alpha is obtained
% This is combined with a local squeezed state or thermal excitation
% Squeezed state parameter is p.sqz, decoherence parameter is p.thermal 
% p.matrix is the transmission, p.phase is the phase-space method = 1,2,3 
% The transmission output can be modified by the factor p.tr(i) 
% Output for +P is a = [alp;bet(;W)], including transmission factors
% Output for Q,W is a = [alp], including transmission factors
% If p.counts > 0, returns a number array, not amplitudes    
% Called by 'ensemble'.
% *Licensed by Peter D. Drummond & Alexander S. Dellios, (2023) 
% - see License.txt, xSDE manual  

%GENERATE PHASE-SPACE SAMPLES FOR ANY ORDERING
 
p   = varargin{end};
w   = varargin{end-1};
sm  = (p.phase-1)/2;                             %%ordering variance
M   = p.modes(1);                                %%abbreviated mode number
r   = p.sqz';                                    %%squeezing parameter
t   = [p.tr,1+0*(1:(M-length(p.tr)))]';          %%transmission 
alp = [p.alpha,0*(1:(M-length(p.alpha)))]';      %%coherent amplitude
n   = sinh(r).^2;                                %%photon number
m   = (1-p.thermal').*cosh(r).*sinh(r) ;         %%initial thermal frac
x   =  sqrt((n+m+sm)/2).*w(1:M,:);               %%x-quadrature noise
y   =  sqrt((n-m+sm)/2).*w(1+M:2*M,:);           %%y-quadrature noise
if p.phase == 1
    bet = (x-1i*y) + conj(alp);                  %%input conjugate noise
end
alp = (x+1i*y) + alp;                            %%input amplitude noise

%ADD OUTPUT LOSSES, TRANSFORMATIONS, EXTRA QUANTUM NOISE

T  = p.matrix(p).*t;                             %%corrected transmission
if isfield(p,'permute') == 1                     %%check for  permutation
    T = T(p.permute,:);                          %%apply random permutation 
end                                              %%end permutation check
alp = T*alp;                                     %%matrix transformation                        
if  sm > 0                                       %%if not normal ordered
  E   = eye(M)-T'*T;                             %%change from unitarity 
  if norm(E) > 1.e-12                            %%if T matrix nonunitary 
    [U,T] = schur(E);                            %%Schur decomposition 
     E1   = U*sqrtm(T)*U';                       %%square root of diagonal
     w1   = (w(1+2*M:3*M,:)+1i*w(1+3*M:4*M,:));  %%additional quantum noise
     w1   = w1*sqrt(sm/2);                       %%quantum noise normalized
     alp  = alp+E1*w1;                           %%output alpha amplitudes
  end                                            %%end if nonunitary
  a   = alp;                                     %%pack Wigner/Q fields
else                                             %%else if normal ordered
  bet = conj(T)*bet;                             %%conjugate transformation
  if p.counts
      a   = alp.*bet;                            %%compute +P number
  else
      a   = [alp;bet];                           %%pack +P fields
  end
end                                              %%end if not normal
end                                              %%end function