function d = SSE(psi,w,p)
% d = SSE(psi,w,p) gives derivatives for a stochastic Schroedinger equation
% The wavefunction is stored as psi(n_1,n_2,n_k...), normalised to unity
% The n-th dissipative operator for mode k is p.L{n}(psi,[k1,..kx],p)
% The n-th amplitude decay rate is p.gamma{n}(k), where k is a mode index.
% If required, k=[k1,..kx] can be a vector index for more than one mode.
% This may be required for losses coupling multiple modes
% The index is passed to the reservoir coupling operator, which is
% EITHER a function p.L{n}(k,psi) for p.sparse = 0, returning a wavefunction
% OR:    a function p.L{n}(k,p)   for p.sparse = 1,  returning an operator
% The non-dissipative Hamiltonian acting on psi is
% EITHER a function p.H(psi,p)    for p.sparse = 0, returning a wavefunction
% OR:    a function p.H(p)        for p.sparse = 1,  returning an operator
% The equation is in the Stratonovich form of stochastic calculus.
%
% Master equn: d rho/dt=-i[H,rho]+sum_n([L_n.rho,L_n^dg]+[L_n,rho.L_n^dg])
% Where: L_n = sum_k(gamma_n(k).L_n(k)
%
% SSE: d psi/dt = -iH. psi + sum_n [(L_n - <L_n>).(z + 2<L_n^dg>)
%                                 + <L_n^dg L_n>-L_n^dg L_n].psi
% where <<dz dz^*>> = 2dt
%
% See: Diosi et al, PRA 58, 1699, (1998) and others, note scaling of gamma.
% For best results, combine this with a projected algorithm like p.MP
%   Input:  wavefunction psi, real noises w, data structure p 
%   Output: derivative d
%   Called by: p.method
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

l = 1;                                           % Initialize noise index
nd = p.nfields;                                  % Store number of modes
if p.sparse                                      % If using sparse methods
  d  = -1i*p.H(p)*psi;                           % Multiply psi by H
else                                             % Else not sparse methods
  d  = -1i*p.H(psi,p);                           % Use H.psi as a function
  zshape = [ones(1,nd),p.ensembles(1)];          % Shape of noise
end                                              % End if sparse methods
for n = p.nl                                     % Loop over loss terms
  gm   = p.gamma{n}(p);                          % Get the n-th loss rate
  szg = size(gm);                                % Get size of loss rate
  for m = 1:numel(gm)                            % Loop on n-th loss vector
    g = gm(m);                                   % Store m-th loss rate
    if g > 0                                     % If loss rate positive
      theta = p.theta{n}(m);                     % Store m-th loss phase
      t2   = theta^2;                            % Store squared phase
      mv = m;                                    % Initialise loss index
      if szg(1)> 1                               % If tensor valued loss
        [mv(1),mv(2)] = ind2sub(szg,m);          % Set tensor indices
      end                                        % End if tensor valued
      if theta == 0                              % If loss phase zero
        z  = sqrt(g)*(w(l,:)+1i*w(l+1,:));       % Set complex noise
        l  = l + 2;                              % Increment noise index
      else                                       % If loss phase nonzero
        z  = sqrt(2*g)*theta*w(l,:);             % Set real noise
        l  = l + 1;                              % Increment noise index
      end                                        % End if loss phase zero
      if p.sparse                                % If sparse method
        L  = p.L{n}(mv,p);                       % Get the loss operator
        X  = L' + t2*L;                          % Get the X operator
        L  = L*psi;                              % Get L*psi
        EX = sum(conj(psi).*(X*psi));            % Get  <X>
        XL = X*L;                                % Get  X*L*psi
      else                                       % Function call method        
        z = reshape(z,zshape);                   % Reshape noise
        L  = p.L{n}(mv,psi);                     % Get L*psi
        X  = p.L{n}(-mv,psi) + t2*L;             % Get X*psi
        XL = p.L{n}(-mv,L) + t2*p.L{n}(mv,L);    % Get  X*L*psi
        EX = sum(conj(psi).*X,1:nd);             % Get  <X>
      end                                        % End if sparse
      d  = d-g*(XL-sum(conj(psi).*XL,1:nd).*psi);% deterministic derivative
      EL = sum(conj(psi).*L,1:nd);               % Get  <X>
      L  = L - EL.*psi;                          % Get  [Delta L]*psi
      d  = d + (2*g*EX+z).*L ;                   % Stochastic derivative  
    end                                          % End if loss rate positive
  end                                            % End loop over loss vector
end                                              % End loop over loss terms
end                                              % End derivative function