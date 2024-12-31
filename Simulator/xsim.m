function [error,data,input,raw] = xsim(varargin)
%  [error,...] = XSIM() solves stochastic differential equations      
%  Solves stochastic differential and partial differential equations.
%  It returns averages as defined by the 'input' cell array. 
%  XSIM also calls parallel loops, estimates errors, writes files.
%  The 'error' is a six dimensional vector:
%   error(1) is the overall RMS (relative) error including differences
%   error(2) is the overall RMS (relative) step error
%   error(3) is the overall RMS (relative) sampling error
%   error(4) is the overall RMS (relative) difference error
%   error(5) is the overall maximum chi-squared error
%   error(6) is the elapsed time
%
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
sx     = max(3,1+p.nchk);                        % sampling error index
error  = zeros(1,6);                             % initialize errors
raw    = cell(seq,p.nchk,p.ensembles(2),p.ensembles(3));
totout = 0;
for s = 1:seq                                    % loop over sequence
  for n = input{s}.outputs                       % loop over outputs
    data{s}{n} = zeros([input{s}.d.data{n},3]);  % initialize output data
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
er = zeros(totout,max(6,p.nchk+4));              % initialize rms errors
o = 0;                                           % initialize total observe
em = zeros(totout,max(6,p.nchk+4));              % initialize max errors
for s = 1:seq                                    % loop over sequence
 p   = input{s};                                 % get input structure
 p.t =  p.origins(1)+(0:p.points{1}(1)-1)*p.dt;  % define time vector
 p.r{1} = p.t;                                   % define first dimension
 for n = p.outputs                               % loop over outputs
  if p.d.data{n}(1) > 0                          % check if data present
    d  = real(data{s}{n});                       % take real part of data
    o = o+1;                                     % increment total count
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
    if p.relerr                                  % If relative errors                      
        mx  = (max(abs(d(:,1))))^2 + 1.e-20;     % max. size for scaling
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
    for nc = 2:p.nchk
      d(:,nc) = abs(d(:,nc) - d(:,1));           % step errors
      er(o,nc) = mean(d(:,nc).^2);
      em(o,nc) = sqrt(max(d(:,nc).^2));
    end
    d(:,2) = sqrt(sum(d(:,2:p.nchk).^2,2));      % RMS sum step errors
    d(:,3) = d(:,sx);
    d      = d(:,1:3);
    esum = (d(:,3).^2+d(:,2).^2)+1.e-20;         % calculate error sum
    ept = (esum > p.tol)&(abs(d(:,1))>p.cutoffs{n});% nonzero errors
    p = p.grid(n,p);                             % get graphics x,o,grid
    p.compares{n} = 0;                           % initialize compares
    
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
      d(:,4:3+ds) = dif; 
      ept = ept&(abs(d(:,4))>p.cutoffs{n});      % Nonzero error points
      mx  = max(abs(d(:,4)).^2); 
      if p.relerr  &&  mx > p.tol                % If relative errors                      
        rl = 'Relative';                         % rel. normalisation
      else                                       % else absolute errors
        mx = 1;
        rl = 'Absolute';                         % abs. normalisation
      end
      p.gpoints{n}(end) = 3+ds;                  % Increase  points
      dif = (d(:,1)-d(:,4)).^2;                  % Get difference square
      em(o,sx+1) = sqrt(max(dif)/mx);            % max mod dif scaled
      er(o,sx+1) = mean(dif)/mx;                 % scale mean dif square
      if ne > 1  
        er(o,sx+2) = sum(ept.*(dif./esum));      % Chi squared
      end
      er(o,sx+3) = sum(ept);                     % Total pts with errors
    end                                          % End check comparisons
    er(o,2:p.nchk) = er(o,2:p.nchk)/mx;
    p.axes{n}  = axes;                           % Compressed axes
    er(o,sx)   = mean(d(:,3).^2)/mx;             % Squared sampling errs
    em(o,sx)   = sqrt(max(d(:,3).^2)/mx);        % Max sampling errors
    er(o,1)    = sum(er(o,2:sx+1))/(1.e-99+sum(er(o,2:sx+1)>0));
    er(o,1:sx+1) = sqrt(er(o,1:sx+1));           % RMS errors
    data{s}{n} = reshape(d,p.gpoints{n});        % Reshape output data
    xpr(1,p,'\n%s errors in output (%d.%d): %s\n',rl, s,n,p.olabels{n});
    for c = 2:p.nchk
      xpr(1,p,'RMS [max] time-step  = %.3g [%.3g]\n',er(o,c),em(o,c));
    end
    xpr(1,p,'RMS [max] sampling   = %.3g [%.3g]\n',er(o,sx),em(o,sx));    
    if ~isequal(p.compare{n},[])
      s1 = sx+1;
      xpr(1,p,'RMS [max] comparison = %.3g [%.3g]\n',er(o,s1),em(o,s1));
      sumpts = mx*1.e-10+er(o,sx+3);
      xpr(1,p,'Chi-square, points  [ratio]  = %.3g, %.3g [%.3g]\n',...
          er(o,sx+2), sumpts, er(o,sx+2)/sumpts);     
    end
    xpr(1,p,'Normalised with max [cutoff] = %.3g [%.3g]\n',sqrt(mx),p.cutoffs{n});
  end                                            % end check if data 
 end                                             % end loop on outputs: n
 if p.dimensions > 1
    p = rmfield(p,{'r','k','x','y','z','kx','ky','kz','t','w'});
 end
 input{s} = p;                                   % store new parameters
end                                              % end loop over sequence
pts = 1.e-99+sum(sum(er(:,2:p.nchk)>p.tol));     % count nonzero datapoints
error(2) = sqrt(sum(sum(er(:,2:p.nchk).^2.))/pts);% average error
error(3) = sqrt(sum(er(:,sx).^2)/(1.e-99+sum(er(:,sx)>p.tol)));
error(4) = sqrt(sum(er(:,sx+1).^2.)/(1.e-99+sum(er(:,sx+1)>p.tol))); 
error(1) = sqrt(sum(error(2:4).^2.)/(1.e-99+sum(error(2:4)>p.tol)));%%store 
error(5:6) = sum(er(:,sx+2:sx+3),1);
error(5)   = error(5)/(1.e-20+error(6));
error(6)   = toc();
pm = functions(p.method);

fprintf('\nRMS Error summary of all outputs\n');
fprintf('\nMethod = %s, Order = %d, Checks = %d',...
          pm.function,p.order,p.checks(1));
fprintf('\nRMS Step+samp+diff error = %.3g\n', error(1));
fprintf('Errors: Step=%.3g Samp=%.3g Diff=%.3g Chisq/k=%.3g\n',error(2:5));
if p.relerr
  fprintf('Errors normalised except for comparisons with zero\n');
else
  fprintf('Using absolute errors\n');
end
fprintf('Computed %d outputs, elapsed time = %.4gs\n',o,error(6));
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