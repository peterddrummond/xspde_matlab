function inputs = xpreferences(inputs)
 % in = XPREFERENCES(inputs) sets default values for the input data.
 % Input:  cell array of parameter structures,'inputs'.
 % Output: cell array of parameter structures, with default values set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Called by xsim
 % Calls qpreferences, phasepreferences, spdpreferences
 % Licensed by Peter D. Drummond, (2023) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sequence = length(inputs);                       % get sequence length
origin = 0;                                      % initial time-origin
for s = 1:sequence                               % loop over sequence 

% PREFERRED VALUES OF PARAMETERS AND FUNCTIONS

  p = inputs{s};                                 % get parameters for s
  p.octave  =     exist('OCTAVE_VERSION', 'builtin') ~= 0;
  p.version =    xprefer(p,'version',0,'xSPDE4.21');
  if ~p.octave
    p.date  =    xprefer(p,'date',0,datetime);   % set the  date in Matlab
  else
    p.date  =    xprefer(p,'date',0,0);          % set the  date in Octave
  end
  p.name =       xprefer(p,'name',0,' ');        % default  name
  p.ensembles =  xprefer(p,'ensembles',3,1);     % default ensembles  
  p.tol       =  xprefer(p,'tol',1,1.e-20);      % tolerance in errors
  p.propfactor = xprefer(p,'propfactor',1,@xpropfactor);
  if s == 1                                      % if first in sequence?
    fprintf('\n%s: %s\n\n',p.version,p.date);    % print version and date
  else                                           % end if first
    p.checks = inputs{1}.checks;                 % ensure matched checks
  end                                            % end if first
  if isfield(p,'quantum') && p.quantum > 0       % if quantum switch?
    p = qpreferences(p);                         % get quantum parameters
  else                                           % or else 
    p.quantum = 0;                               % set p.quantum = 0
  end                                            % end if quantum?
  if isfield(p,'phase') && p.phase > 0           % if phase-space?
    p = phasepreferences(p);                     % phase-space parameters
  end                                            % end if phase-space
  p.dimensions = xprefer(p,'dimensions',1,1);    % partial diff. switch
  nd = max(1,p.dimensions);                      % space-time dimensions
  p.fields =     xprefer(p,'fields',0,{1});      % integrated fields
  p.fieldsb =    xprefer(p,'fieldsb',0,[]);      % backfields
  p.auxfields =  xprefer(p,'auxfields',0,[]);    % auxiliary fields  
  p.auxcells   = length(p.auxfields);            % number of auxil. cells
  p.fieldcells = length(p.fields);               % number of forward cells
  p.fbcells    = length(p.fieldsb)+p.fieldcells; % number of total cells
  p.fields     = [p.fields,p.fieldsb];           % total fields
  p.totcells   = p.fbcells+p.auxcells;           % total cells
  if p.dimensions == 0                           % if no space-time  
    p.points   = xprefer(p,'points',1,{1});      % only one time point
    p.checks   = xprefer(p,'checks',0,0);        % no error checking
  else                                           % else has dimensions 
    p.points   = xprefer(p,'points',1,{[51,35*ones(1,nd-1)]}); % default
    p.checks   = xprefer(p,'checks',0,[1,zeros(1,nd-1)]);% error checks
  end                                            % end if no space-time
  p.iterations = xprefer(p,'iterations',1,4);    % number of iterations
  p.iterfb     = xprefer(p,'iterfb',1,1);        % forward-backward iters
  p.breed  =     xprefer(p,'breed',1,@Breed);    % default breeding method
  p.path   =     xprefer(p,'path',1,@xpath);     % default SDE path
  p.pathfb   =   xprefer(p,'pathfb',1,@xpathfb); % default FBSDE path
  p.deriv      = xprefer(p,'deriv',0,{@(a,~,~) 0*a});
  p.noises     = xprefer(p,'noises',0,{[]});     % gauss noises
  p.unoises    = xprefer(p,'unoises',0,{[]});    % uniform noises
  p.knoises    = xprefer(p,'knoises',0,{[]});    % filtered noises
  if isequal({p.noises,p.knoises,p.unoises},{{[]},{[]},{[]}}) % if null
    p.noises   = p.fields;                       % use default noises
  end                                            % end if null  
  p.inrandoms  = xprefer(p,'inrandoms',0,{[]});  % gauss randoms                                          
  p.urandoms   = xprefer(p,'urandoms',0,{[]});   % uniform random
  p.krandoms   = xprefer(p,'krandoms',0,{[]});   % filtered randoms
  if isequal({p.inrandoms,p.krandoms,p.urandoms},{{[]},{[]},{[]}})%if 0
    p.inrandoms = p.noises;                      % use default inrandoms
  end                                            % end if null
  p.vc         = xprefer(p,'vc',0,{[]});         % polynomial vector coeff.
  p.orderpol   = xprefer(p,'orderpol',0,{[]});   % polynomial order
  p.qc         = xprefer(p,'qc',0,{[]});         % quadratic matrix coeff.
  p.fields     = [p.fields,p.auxfields];         % total fields + auxfields
  p.ranges     = xprefer(p,'ranges',1,10*ones(1,nd));
  p.steps =      xprefer(p,'steps',nd,ones(1,nd));
  p.origins =    xprefer(p,'origins',1,[origin,-p.ranges(2:nd)/2]);   
  p.jump   =     xprefer(p,'jump',p.fieldcells,{0});
  p.rawdata =    xprefer(p,'rawdata',0,0);       % switch for raw data
  p.linear   =   xprefer(p,'linear',p.fbcells,{0});% linear prop.
  p.adapt   =    xprefer(p,'adapt',1,1);
  p.seed =       xprefer(p,'seed',1,0);
  p.rng =        xprefer(p,'rng',0,'Threefry');
  p.file =       xprefer(p,'file',0,' ');
  p.verbose =    xprefer(p,'verbose',1,0);   
  p.initial =    xprefer(p,'initial',0, {@(~,~) 0;});
  p.define  =    xprefer(p,'define',0,{@(varargin) []});
  p.transfer =   xprefer(p,'transfer',0,{@(a,~,~) a});
  p.thresholdw = xprefer(p,'thresholdw',1,0);
  p.breedw  =    xprefer(p,'breedw',1,0);        % default breeding
  p.relerr   =   xprefer(p,'relerr',1,1);        % switch for rel. errors
  p.counts   =   xprefer(p,'counts',1,0);        % number of jumps 
  p.spectrum =   xprefer(p,'spectrum',1,0);
  p.grid =       xprefer(p,'grid',0,@xgrid);     % handle for grid
  if isequal(prod(p.ensembles),1)  
      p.method = xprefer(p,'method',1,@RK4);     % preferred deterministic
    else 
      p.method = xprefer(p,'method',1,@MP);      % preferred stochastic
  end
  if p.iterfb > 0
      p.method = xprefer(p,'method',1,@MPfb);  % preferred fb method
  end
  p.firstfb =  xprefer(p,'firstfb',p.fbcells,{@(p) 0});     % preferred first fb
  p.prop     =   xprefer(p,'prop',1,@xprop);     % preferred propagator
  
                                        
% CALCULATE CONSTANTS THAT DEPEND ON THE CELL INDEX

  p.nfields      = 0;
  for c = 1:p.fieldcells                         % index over field cells
    p.nfields = max(length(p.fields{c}),p.nfields);
  end
  p.ind(1:p.dimensions+p.nfields+1) = {':'};     % store index cells
  p.noisegen =   xprefer(p,'noisegen',1,@xnoise);
  p.randomgen =  xprefer(p,'randomgen',1,@xrandom);
  p.nfilter =    xprefer(p,'nfilter',p.fieldcells,{@xnfilter});
  p.rfilter =    xprefer(p,'rfilter',p.fieldcells,{@xrfilter});
  dx = cell(1,p.totcells);   a = dx; space = dx;
  ens = p.ensembles(1);
  p.propagator = num2cell(ones(p.fieldcells,nd));
  p.da = num2cell(zeros(p.fieldcells,nd));
  mc = max([p.totcells,numel(p.noises),numel(p.knoises),numel(p.unoises)...
       ,numel(p.inrandoms),numel(p.krandoms),numel(p.urandoms)]);

  % CALCULATE CELL PARAMETERS
  
  p.points(end+1:mc) = p.points(end);
  for c = 1:mc                                   % index over all cells
    p.points{c}(1)    = p.points{1}(1);
    p.points{c}       = [p.points{c},ones(1,nd-length(p.points{c}))];  
    p.points{c}(2:nd) = 1+(p.points{c}(2:nd)-1).*p.steps(2:nd);
    dx{c} = p.ranges./max(1,p.points{c}-1);      % n-th step-size in x
    dx{c} = (p.points{c} == 1)+(p.points{c} > 1).*dx{c};
    p.dv{c}  = prod(dx{c}(2:nd));                % lattice cell volume 
    p.noisefactor{c}  = 1/p.dv{c};
    p.knoisefactor{c} = prod(p.points{c}(2:nd))/p.dv{c};
    space{c}   = p.points{c}(2:nd);
    if c <= p.totcells
      p.d.ca{c}  = [p.fields{c},space{c},ens];
      p.d.a1{c}  = [p.fields{c},1,space{c},ens];
      p.d.ft{c}  = [p.fields{c},p.points{c}*p.steps(1),ens];
      p.d.raw{c} = [p.fields{c},p.points{c},ens]*p.rawdata;
      p.d.rc{c}  = [prod(p.fields{c}),p.points{c},ens];
      a{c}       = zeros(p.d.ca{c});          % Initial data 
    end
  end
  p.d.rand     = cappend(p.inrandoms,space,ens);
  p.d.krand    = cappend(p.krandoms,space,ens);
  p.d.urand    = cappend(p.urandoms,space,ens);
  p.d.noise    = cappend(p.noises,space,ens);
  p.d.knoise   = cappend(p.knoises,space,ens);
  p.d.unoise   = cappend(p.unoises,space,ens);
  p.noisecells = [numel(p.d.noise),numel(p.d.knoise),numel(p.d.unoise)];
  p.randcells  = [numel(p.d.rand),numel(p.d.krand),numel(p.d.urand)];
  p.d.r        = [ones(1,p.nfields),p.points{1}(2:nd)]; 
  p.indext     = 0;
  p.dx         = dx{1};
  p.dv         = p.dv{1};
  p.dk = 2.0*pi./(p.points{1}.*p.dx);            % n-th step-size in k
  p.dk(1) = p.dk(1)/p.steps(1);                  % 1st step-size in k
  p.kfact =     p.dx/sqrt(2*pi);                 % fft (k) normalization
  p.nspace =     prod(p.points{1}(2:end));       % count space points
  if nd > 1 
    p = spdpreferences(p);                       % Partial differential
  else
    p.setboundaries = 0;
  end
  p.dt      = p.ranges(1)/(max(1,p.points{1}(1)-1));
  p.nscale  = p.steps(1)*(1+p.checks(1))/p.dt;   % noise scale factor
  if isequal (p.linear,{[]}) && nd == 1
      p.prop    =   @(a,~,~) a;
  end
  nc = 1;
  p.checkd(1) = 0;
  p.spacechecks = sum(p.checks(2:nd));
  p.nchk = 1+sum(p.checks);
  p.fbsw = abs(p.iterfb)>1;
  for n = 1:nd
    if p.checks(n)
      nc = nc+1;
      p.checkd(nc) = n;
      p.dxr(n) = p.dx(n)/(2*p.steps(n)); 
    else
      p.dxr(n) = p.dx(n)/(p.steps(n));
    end
  end
  
% CALCULATE GENERAL CONSTANTS, 

  p.d.a       =  p.d.ca{1};                      % Dimension of first cell
  origin      =  p.origins(1) + p.ranges(1);     % Origin of next sequence
  p.errors    =  3;                              % Maximum error field
  p.dw        =  2*pi/(p.ranges(1)+p.dt);
  p.fsc       =  sqrt(2*pi)/p.dw;                % Fourier scale factor
  p0          =  (0:p.points{1}(1)-1);
  p.plotpts   =  (p0-floor((p.points{1}(1))/2));
  p.w         =  p.plotpts*p.dw;  
  nspace      =  prod(space{1});
  time        =  p.points{1}(1);
  p.xc{1}     =  p.origins(1) + p.dx(1)*p0;      % 1st x-coords
  p.kc{1}     =  p.w;                            % 1st k-coords
  if nd > 1 
      p = p.grid(0,p);                           % get grid-points
  end
  p.dtr = 0; 
  p.t = p.origins(1);
  if isequal (p.method, @Euler) || isequal (p.method, @Implicit)...
      || isequal (p.method, @RK2)    
      p.ipsteps = xprefer(p,'ipsteps',1,1);
  else 
      p.ipsteps = xprefer(p,'ipsteps',1,2);
  end
  p.order = xprefer(p,'order',1,0);
 
  
% GET PARAMETERS FOR THE OUTPUT AVERAGES

  if ~isfield(p,'observe')
    p.observe =  {@xobserve};      % set default observe
  else
    p.observe =  xmcell(p.observe);              % make observe a cell
  end
  p.expect    =  xprefer(p,'expect',0,{[]});     % set default expect
  for  e = 1:length(p.expect)
    if ~isequal(p.expect{e},[]) && p.thresholdw  % make weighted expect
      p.observe{e} = @(a,om,p) p.expect{e}(a,p).*exp(om)./mean(exp(om));                   
    end
  end
  if ~isfield(p,'averages')
      p.averages = 1:length(p.observe);
  end
  p.scatters  =  xprefer(p,'scatters',max(p.averages),{0});
  p.binranges =  xprefer(p,'binranges',max(p.averages),{{}});% bin range
  p.transforms = xprefer(p,'transforms',max(p.averages),{zeros(1,nd+1)});
  p.mincount  =  xprefer(p,'mincount',1,10);
  p.transforms(end+1) = {0};  
  O = cell(1,max(p.averages));
  p.d.av = cell(1,max(p.averages));
  for n = p.averages  
    p.spectrum =   p.spectrum || p.transforms{n}(1) ;
    for j = 1:nd
      if ~p.transforms{n}(j)
        p.xk{n}{j}  = p.origins(j)+(0:p.points{1}(j)-1)*p.dx(j);
      else
        start = -p.points{1}(j)/2;
        p.xk{n}{j}  = (start:start+p.points{1}(j)-1)*p.dk(j);
     end
    end
  end
  
  fprintf('Sequence %d, %s\n',s,p.name);         % print the sequence id
  fprintf('Computing %d average(s)\n',numel(p.averages));% print averages count
 
% CHECK OBSERVE DIMENSIONS FROM TEST INPUTS

  for n = p.averages                             % loop over averages
    O{n} = [];                                   % set observe to void
    if ~isempty(p.observe{n})                    % if observe defined
      p.nobserve = n;
      O{n} = p.observe{n}(a{:},p);               % return observe output
      sz = size(O{n});                           % size of observe output
      szl = sz(end);                             % get last observe size
      if szl > ens && rem(szl,ens) == 0          % if last size > ens
          sz(end)=szl/ens;                       % factor: stop mismatch
          sz(end+1)=ens;                         %#ok<AGROW>
      end                                        % end check last size
      switch p.dimensions
          case 0
              lines = 1;                         % get number of lines
              spaceobs = sz(1:end-1);            % get space points
              ens1 = sz(end);                    % get ensemble points
              p.xk{n}{1} = 1;
              for m = 1:numel(spaceobs)
                   p.xk{n}{1+m}  = 1:spaceobs(m);
              end
          case 1
              lines = sz(1:end-1);               % get number of lines
              spaceobs = [];                     % get space points
              ens1 = sz(end);
          otherwise                              % if space dimension
              lines = sz(1);                     % get number of lines
              if prod(sz(2:end-1)) == nspace
                  spaceobs = space{1};
                  ens1 = sz(end);
              else
                  spaceobs = sz(2:end-1);        % get space points
                  ens1 = sz(end);
              end
      end                                        % end if space
      p.d.obs{n} = [prod(lines),1,prod(spaceobs),ens1];
      spacebins = spaceobs;

% CALCULATE SCATTER AND PROBABILITY SIZES

      p.scatters{n} = min(p.scatters{n},sz(end));
      if p.scatters{n} > 1                       % check for scatters 
        lines = lines*p.scatters{n};             % modify lines
        if p.ensembles(3) > 1                    % If parallel
          fprintf ('\nWARNING: scatters{%d}+parallels is incompatible',n);
          fprintf ('\nResetting parallels to 1\n');
          p.ensembles(3) = 1;
        end                                      % end if verbose printing
      end                                        % end check for scatters
      p.d.bins{n} = [];                          % Initial dimension
      binranges  = p.binranges{n};               % n-th bin ranges
      p.bins{n} = 0;                             % n-th bin numbers
      if ~isempty(binranges)                     % If n-th bins set         
        p.bins{n} = min(length(binranges),lines);
        dbin = zeros(1,p.bins{n});
        abin = 1;                                % Initial area of bin;
        m = 0;
        for m1 = 1:p.bins{n}                     % Loop on probability axes 
          d = length(binranges{m1})-1;
          if d > 0
             m = m+1;
             dbin(m) = d;                        % Dimension of n-th bin
             p.oc{n}{m} = (binranges{m1}(2:d+1)+binranges{m1}(1:d))/2;
             p.do{n}{m} = (binranges{m1}(2:d+1)-binranges{m1}(1:d));
             abin = abin*p.do{n}{m}(1);          % Area of n-th bin          
             p.xk{n}{nd+m}  = p.oc{n}{m};        % Make a new axis vector
          end
        end                                      % End loop over length
        p.bins{n} = m;                           % n-th bin numbers
        if m > 0
          if p.verbose > 0                       % If verbose printing
            fprintf ('\n%d#%d is a %d D probability',s,n,p.bins{n});
          end                                    % end if verbose printing
          p.d.bins{n} = dbin;                    % Dimensions of bin
          p.a.bins{n} = abin;                    % Area of bin
          spacebins = [spaceobs,dbin];           % Dimension of space+bins
          lines = 1;                             % Lines = 1: probability
        end
      end                                        % End if probability bins
      p.d.av{n} = [prod(lines),time,spacebins];  % Dimension of averages
      O{n} = zeros(p.d.av{n});                   % Observe test array cell
    end                                          % End if observe present
  end                                            % End loop over averages
  
% CALCULATE OUTPUT PARAMETERS AND SIZES
    
  p.output =   xprefer(p,'output',0,{[]});
  for n = p.averages
    if n > length(p.output) || isequal(p.output{n},[]) 
      p.output{n} = @(o,p) o{n};
    end 
  end
  i = 1;
  for n = 1:length(p.output)
      if ~isequal(p.output{n},[]) 
        p.outputs(i) = n;
        i=i+1;
      end
  end
  p.olabels =    xprefer(p,'olabels',max(p.outputs),{''});
  axes{1} = num2cell(zeros(1,nd+10));  % initialize default axes
  fprintf('\nCalculating %d output(s)\n',numel(p.outputs));
  p.compare =     xprefer(p,'compare',max(p.outputs),{[]});
  p.ftransforms = xprefer(p,'ftransforms',max(p.outputs),p.transforms);
  p.binranges =   xprefer(p,'binranges',max(p.outputs),{{}});
  p.cutoff  =   xprefer(p,'cutoff',0,-1.e20);
  p.cutoffs  =   xprefer(p,'cutoffs',max(p.outputs),{p.cutoff});
  for n = p.outputs                            % get output sizes
    if n>max(p.averages)
      p.xk{n} = p.xk{max(p.averages)};           % get unset coordinates
      p.scatters{n} = 0;
      p.bins{n} = 0;
    end
    if isequal (p.olabels{n}, [])
        p.olabels{n} = ' ';
    end
    p.spectrum = p.spectrum || p.ftransforms{n}(1);
    p.gpoints{n} = size(p.output{n}(O,p));
    p.d.data{n}  = prod(p.gpoints{n});
    p.gpoints{n}(1+end) = max(3,1+p.nchk+p.fbsw);
    p.xk{n}{1}  = p.xc{1};
    p.plotpts = ((0:p.points{1}(1)-1)-floor((p.points{1}(1))/2));
    if ~isfield(p,'axes') || numel(p.axes) < n 
      for j= 2:length(p.steps)
        if p.steps(j) > 1 
          if ~p.ftransforms{n}(j)
            axes{1}{j} = 1:p.steps(j):p.points{1}(j);
          end
        end
      end
    end
    for j=1:nd                         %loop on dimensions
      if p.ftransforms{n}(j)
        points = p.points{1}(j)-1;
        if ~isfield(p,'axes') || numel(p.axes) < n 
          if j > 1 && j <= length(p.steps) && p.steps(j)>1
            points = points/p.steps(j);
            mid = (p.points{1}(j)-points)/2;
            axes{1}{j} = mid:mid+points;
          end
        end
        start = -points/2;
        p.xk{n}{j}   = (start:start+points)*p.dk(j);
      end
    end                                          %end loop on dimensions
  end                                            %end loop on outputs
  axes{1} = num2cell(zeros(1,nd+10));
  p.axes = xprefer(p,'axes',max(p.outputs),axes);
  inputs{s} = p;
  if p.verbose > 1
    display(p);
    display(p.d);
  end
end                                              % end loop on sequence
end                                              % end xpreferences

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END XPREFERENCES FUNCTION

function kn = xnfilter(w,~)                      % Default noise filters
    kn = w;                                      % Default noise filter
end

function kr = xrfilter(w,~)                      % Default random filters
    kr = w;                                      % Default input filter 
end

function o = xobserve(varargin)                 % Default observe
% O = XOBSERVE(varargin) is the default xSPDE observe function
 % Input:  field cells, parameters
 % Output: Default observe: all field components as separate lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Called by xpreferences
 % Licensed by Peter D. Drummond, (2024) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    o   = varargin{1};                          % first input cell
    p   = varargin{end};                        % parameters
    sz  = size(o);
    elin = numel(sz);
    el  = numel(p.fields{1});
    if elin <= el
      szo = [prod(sz(1:elin)),1];
    else
      szo = [prod(sz(1:el)),sz(el+1:end)];
    end
    o   = reshape(o,szo);
end

function cout = cappend(cin,cvec,vec)
 % COUT = CAPPEND(CIN,CVEC,VEC) appends vectors to cells of vectors
 % Input:  cells cin, cells cvec with vectors, and a third vector, vec.
 % Output: cells cout{i} = [cin{i} ,cvec{i}, vec] of appended vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Called by xpreferences
 % Licensed by Peter D. Drummond, (2024) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cin = xmcell(cin);                               % Makes cell if not null
cout = [];
if iscell(cin) && ~isequal(cin,{[]})
  cinl  = length(cin);
  cvecl = length(cvec);
  cout = cell(1,cinl);
  for i = 1:cinl
    cout{i} = [cin{i},cvec{min(i,cvecl)},vec];    % Extend i-th vector
  end
end
end