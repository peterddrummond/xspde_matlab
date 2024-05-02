function p = qpreferences(p)
%   p = QPREFERENCES(p) sets default values for quantum input data.
%   Input:  original parameter structure,'p'.
%   Output: parameter structure, with quantum default values set.
%   Called by xpreferences
%   Licensed by Peter D. Drummond, (2024) - see License.txt, xSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.nmax    = xprefer(p,'nmax',0,2);               % default mode dimension
p.modes   = xprefer(p,'modes',0,numel(p.nmax));  % default mode number
p.sparse  = xprefer(p,'sparse',0,0);             % default sparse switch
p.jump    = xprefer(p,'jump',0,0);               % default jump switch
p.gamma   = xprefer(p,'gamma',1,{@(p) 0});       % default loss function
p.L       = xprefer(p,'L',1,{[]});               % default loss operator
p.theta   = xprefer(p,'theta',1,{0});            % default alpha
p.measure = xprefer(p,'measure',1,0);            % measured terms in SME
nm = min(p.modes,length(p.nmax));                % test length of nmax
p.nmax(1+nm:p.modes) = p.nmax(nm);               % set nmax on all modes
p.fields = prod(p.nmax);                         % sparse field dimension
p.nl  = 1:min(numel(p.L),numel(p.gamma));        % decay channels simulated
p.t = 1.e-15;                                    % set time to small value
l = 0;                                           % initialize random count
p.theta{max(p.nl)+1}=[];                         % ensure alpha cell exists
for n = p.nl                                     % loop over decay channels
  g = p.gamma{n}(p);                             % get decay rate
  for m = 1:numel(g)                             % loop over decay rates
    if m > numel(p.theta{n})                     % test if alpha is defined
      p.theta{n}(m) = 0;                         % set alpha to default
    end                                          % end test if alpha
    if g(m) ~= 0.0 && (p.quantum < 2 || n <= p.measure)%%test nonzero decay   
      l = l+1;                                   % increment random count
      if p.theta{n}(m) == 0 && p.jump == 0       % test for complex noise
        l = l+1;                                 % increment random count
      end                                        % end test for complex
    end                                          % end test if nonzero 
  end                                            % end loop over decays
end                                              % end loop over channels
if p.jump == 1
    p.unoises = xprefer(p,'unoises',0,{l});      % default uniform noise
    p.noises = xprefer(p,'noises',0,[]);         % default gaussian noise 
else
    p.unoises = xprefer(p,'unoises',0,[]);       % default uniform noise
    p.noises = xprefer(p,'noises',0,{l});        % default gaussian noise 
end
if p.quantum == 2                                % check if master equation
    p.sparse = 1;                                % only sparse is allowed
    p.deriv  =   xprefer(p,'deriv',0,@ME);       % set derivative
    p.fields(2) = p.fields(1);                   % make a square matrix
    p.ensembles(1) = 1;                          % no local ensembles
elseif  p.quantum == 1                           % check if SSE
    if p.jump == 1                               % check if jump
       p.jmethod = xprefer(p,'jmethod',0,@JSEB); % set the jump method
       p.deriv  =  xprefer(p,'deriv',0,@JSEA);
    else                                         % else continuous case
       p.deriv  =   xprefer(p,'deriv',0,@SSE);
    end
end
if p.sparse                                      % check if sparse
    p.fields = {p.fields};
    p.vac = mkvacs(p);
    p.H   = xprefer(p,'H',0, @(p) sparse(p.fields{1}(1),p.fields{1}(1)));
    p.L   = xprefer(p,'L',0, {[]});
else                                             % not sparse
    p.fields = {p.nmax};                         % full field dimension
    p.H = xprefer(p,'H',0,@(a,p) zeros([p.fields,p.ensembles(1)]));
    p.L = xprefer(p,'L',0, {[]});
end  
p.observe  =  xprefer(p,'observe',0,{[]});
p.expect  =   xprefer(p,'expect',0,{[]});
for  e = 1:length(p.expect)
  if ~isequal(p.expect{e},[])
  switch p.sparse+2*(p.quantum - 1)
    case 0
      p.observe{e} = @(psi,p) reshape(sum(conj(psi).*...
      p.expect{e}(psi,p),1:p.modes),[1,size(psi,p.nfields+1)]);
    case 1
      p.observe{e} = @(psi,p) sum(conj(psi).*(p.expect{e}(p)*psi));
    case 2
      error('Only sparse master equations currently suppported');
    case 3
      p.observe{e} = @(rho,p) reshape(...
      trace(p.expect{e}(p)*rho),[1,1]);
      if p.jump
        error('Stochastic master equations not currently suppported');
      end
  end
  end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XQPREFERENCES FUNCTION 