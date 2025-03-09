function [error,data,input,raw] = xsim(varargin)
%  [error,...] = XSIM() solves stochastic differential equations      
%  Solves stochastic differential and partial differential equations.
%  It returns averages as defined by the 'input' cell array. 
%  XSIM also calls parallel loops, estimates errors, writes files.
%  Error summaries include n > 0 step/iteration errors, where:
%  
%   error(1) is the overall RMS (relative) error including differences
%   error(2:1+n) is the RMS (relative) step/iteration error
%   error(2+n) is the overall RMS (relative) sampling error
%   error(3+n) is the overall RMS (relative) difference error
%   error(4+n) is the overall maximum chi-squared error
%   error(5+n) is the overall run-time
%
%  Note that 'error' is a vector of at least 6 dimensions
%  
%  Error(2) is zero if there is no error checking or iterations
%  Error(1+n) stores forward-backward iteration errors, where n > 1
%  The 'data' cell array   gives the averages, errors, and comparisons.
%  The 'input' cell array  gives input data including all default values.
%  The 'raw'  cell array   gives raw trajectories if required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Called by: xspde
%  Calls: xmakecell, xpreferences, xensemble, addcell, xwrite
%  Licensed by Peter D. Drummond, (2024) - see License.txt, xSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic();                                           % start timer
if iscell(varargin{1})
    varargin = varargin{1};
end
input  = xpreferences(varargin);                 % set default parametera
p      = input{1};                               % get first parameters
seq    = length(input);                          % store sequence length
data   = cell(seq);                              % initialize data cells
ne     = p.ensembles(2)*p.ensembles(3);          % total higher ensembles
ens    = prod(p.ensembles);                      % initialize ensembles
n1     = ne - 1;
nc1    = p.nchk+p.fbsw;
sx     = max(3,1+nc1);                           % sampling error index
error  = zeros(1,6);                             % initialize errors
raw    = cell(seq,p.nchk,p.ensembles(2),p.ensembles(3));
totout = 0;
for s = 1:seq                                    % loop over sequence
  for n = input{s}.outputs                       % loop over outputs
    data{s}{n} = zeros([input{s}.d.data{n},sx]); % initialize output data
  end                                            % end loop over outputs
  totout = totout+max(input{s}.outputs);
end                                              % end loop over sequence
if p.ensembles(3) > 1                            % if parallel simulations
  parfor npar = 1:p.ensembles(3)                 % loop on parallel threads
    tic();                                       % set timer
    [dp,raw1] = xensemble(npar,input);           % call xensemble
    raw(:,:,:,npar) = raw1;                      % accumulate raw data
    data = addcell(data,dp);                     % accumulate averages
  end                                            % end parallels loop
  raw = reshape(raw,[seq,p.nchk,ne]);            % reshape raw 
else                                             % no parallels specified
  [data,raw] = xensemble(1,input);               % call xensemble
end                                              % end if parallels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Calculate final mean and variance for sampling errors            % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er = zeros(totout,max(6,nc1+4));                 % initialize rms errors
o = 0;                                           % initialize total observe
em = zeros(totout,max(6,nc1+4));                 % initialize max errors
for s = 1:seq                                    % loop over sequence
 p   = input{s};                                 % get input structure
 p.t =  p.origins(1)+(0:p.points{1}(1)-1)*p.dt;  % define time vector
 p.r{1} = p.t;                                   % define first dimension
 for n = p.outputs                               % loop over outputs
  if p.d.data{n}(1) > 0                          % check if data present
    d  = real(data{s}{n});                       % take real part of data
    o = o+1;                                     % increment output count
    d = reshape(d,p.gpoints{n});                 % output: full indexing
    m = min(p.dimensions,length(p.gpoints{n})-2);% Fourier dimensions
    for nd=2:m                                   % loop over dimensions
      if p.ftransforms{n}(nd)                    % check Fourier transforms
         d = fftshift(d,1+nd);                   % shift Fourier origin
      end                                        % end check Fourier
    end                                          % end loop over dimensions
    [d,p.xk{n},axes] = xcompress(n,d,p);         % Compress data to axes
    p.gpoints{n} = size(d);
    p.d.data{n} = prod(p.gpoints{n}(1:end-1));
    gptlen = length(p.gpoints{n});
    d = reshape(d,[p.d.data{n},sx]);
    mx  = (max(abs(d(:,1))))^2;                  % max. size for scaling
    if p.relerr && mx > p.tol                    % If relative errors                      
        rl = 'Relative';                         % rel. normalisation
    else                                         % else absolute errors
        mx = 1;                                  % then don't scale
        rl = 'Absolute';                         % abs. normalisation
    end                                          % end if relative 
    if ne>1
      d(:,sx)=abs(sqrt((d(:,sx)-d(:,1).^2)/n1)); %  sampling error
    else
      d(:,sx) = 0*d(:,sx);
    end
    if p.bins{n} > 0                             % If probabilities
       p.scale{n} = p.a.bins{n}*ens;             % scale variance
       d(:,sx) = sqrt(d(:,1)/p.scale{n});        % Poissonian sd
       p.cutoffs{n} = max(p.cutoffs{n},p.mincount/p.scale{n});
    end                                          % end if probabilities
    for nc = 2:(p.nchk+p.fbsw)
      d(:,nc)  = abs(d(:,nc) - d(:,1));          % step or iteration errors
      er(o,nc) = mean(d(:,nc).^2);
      em(o,nc) = sqrt(max(d(:,nc).^2));
    end
    esum = sum(d(:,2:sx).^2,2)+1.e-20;           % calculate error sum
    ept = (esum > p.tol)&(abs(d(:,1))>p.cutoffs{n});% nonzero errors
    p = p.grid(n,p);                             % get graphics x,o,grid
    p.compares{n} = 0;                           % initialize compares
    ds = 0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Calculate comparison data points if present                      % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if ~isequal(p.compare{n},[])                 % Check for comparisons
      p.noutput = n;
      dif = p.compare{n}(p);                     % Compute comparison
      dif = dif+zeros([p.gpoints{n}(1:(end-1)),1]);  % Pad zero to size
      [dif,~,~] = xcompress(n,dif,p);            % Compress data to axes
      ds = size(dif,gptlen);
      dif = reshape(dif,[p.d.data{n},1:ds]);     % get the comparisons
      d(:,sx+1:sx+ds) = dif; 
      ept = ept&(abs(d(:,sx+1))>p.cutoffs{n});   % Nonzero error points
      mx  = max(abs(d(:,sx+1)).^2); 
      if p.relerr  &&  mx > p.tol                % If relative errors                      
        rl = 'Relative';                         % rel. normalisation
      else                                       % else absolute errors
        mx = 1;
        rl = 'Absolute';                         % abs. normalisation
      end
      p.gpoints{n}(end) = sx+ds;                 % Increase  points
      dif = (d(:,1)-d(:,sx+1)).^2;               % Get difference square
      em(o,sx+1) = sqrt(max(dif)/mx);            % max mod dif scaled
      er(o,sx+1) = mean(dif)/mx;                 % scale mean dif square
      if ne > 1  
        er(o,sx+2) = sum(ept.*(dif./esum));      % Chi squared
      end
      er(o,sx+3) = sum(ept);                     % Total pts with errors
    end                                          % End check comparisons
    er(o,2:p.nchk) = er(o,2:p.nchk)/mx;
    p.axes{n}  = axes;                           % Compressed axes
    er(o,sx)   = mean(d(:,sx).^2)/mx;             % Squared sampling errs
    em(o,sx)   = sqrt(max(d(:,sx).^2)/mx);        % Max sampling errors
    er(o,1)    = sum(er(o,2:sx+1))/(1.e-99+sum(er(o,2:sx+1)>0));
    er(o,1:sx+1) = sqrt(er(o,1:sx+1));           % RMS errors
    d(:,2) = sqrt(sum(d(:,2:nc1).^2,2));         % RMS sum step errors
    d(:,3) = d(:,sx);
    p.gpoints{n}(end) = 3;                       % Increase  points
    if ~isequal(p.compare{n},[])
      d(:,4:3+ds) = d(:,sx+1:sx+ds);
      p.gpoints{n}(end) = 3+ds;
    end
    data{s}{n}=reshape(d(:,1:3+ds),p.gpoints{n});% Reshape output data
    xpr(1,p,'\n%s errors in output (%d.%d): %s\n',rl, s,n,p.olabels{n});
    for c = 2:p.nchk
      e = c-1;
      xpr(1,p,'RMS[max] step (D%d)  = %.3g [%.3g]\n',e,er(o,c),em(o,c));
    end
    if p.fbsw
      c = p.nchk+1;
      xpr(1,p,'RMS[max] iteration  = %.3g [%.3g]\n',er(o,c),em(o,c));
    end
    xpr(1,p,'RMS[max] sampling   = %.3g [%.3g]\n',er(o,sx),em(o,sx));    
    if ~isequal(p.compare{n},[])
      s1 = sx+1;
      xpr(1,p,'RMS [max] comparison = %.3g [%.3g]\n',er(o,s1),em(o,s1));
      sumpts = mx*1.e-10+er(o,sx+3);
      xpr(1,p,'Chi-square, points  [ratio]  = %.3g, %.3g [%.3g]\n',...
          er(o,sx+2), sumpts, er(o,sx+2)/sumpts);     
    end
    xpr(1,p,'Normalised by %.3g [cutoff = %.3g]\n',sqrt(mx),p.cutoffs{n});
  end                                            % end check if data 
 end                                             % end loop on outputs: n
 p = rmfield(p,{'r','k','x','y','z','kx','ky','kz','t','w','dv'});
 input{s} = p;                                   % store new parameters
end                                              % end loop over sequence
for c = 2:sx+1
  pts = 1.e-99+sum(er(:,c)>p.tol);               % nonzero datapoints
  error(c) = sqrt(sum(er(:,c).^2.)/pts);         % average error
end
error(1) = sqrt(sum(error(2:sx+1).^2.)/(1.e-99+sum(error(2:sx+1)>p.tol)));
error(sx+2:sx+3) = sum(er(:,sx+2:sx+3),1);
error(sx+2)   = error(sx+2)/(1.e-20+error(sx+3));
error(sx+3)   = toc();
pm = functions(p.method);
fprintf('\nRMS Error summary of all outputs\n');
fprintf('\nMethod = %s, Order = %d, Checks = %d',...
          pm.function,p.order,p.checks(1));
fprintf('\nRMS average of errors + differences = %.3g\n', error(1));
for c = 2:p.nchk
  fprintf('Step errors (Dimension % d) = %.3g\n',c-1, error(c));
end
if p.fbsw
  fprintf('Iteration errors = %.3g\n',error(nc1));
end
fprintf('Sampling errors=%.3g Diffs =%.3g Chisq/k=%.3g\n',error(sx:sx+2));
if p.relerr
  fprintf('Errors normalised except for comparisons with zero\n');
else
  fprintf('Using absolute errors\n');
end
fprintf('Computed %d outputs, elapsed time = %.4gs\n',o,error(sx+3));
filename = input{1}.file;
if ~strcmp(filename,' ')
  save (filename,'error','data','input','raw'); % save  Matlab file
end                                             % end if file not blank
if isscalar(input)                              % if  input length is 1
    input = input{1};                           % return struct, not cell
end                                             % end if input length is 1
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data1 = addcell(data1,data2)
%  data1 = ADDCELL() adds two cell array components             
%  Called by sim
%  Licensed by Peter D. Drummond, (2023) - see License.txt, E manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for seq = 1:length(data1)
  for n = 1:length(data1{seq})
    data1{seq}{n} = data1{seq}{n}+data2{seq}{n};
  end
end
end