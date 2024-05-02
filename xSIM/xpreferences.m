function inputs = xpreferences(inputs)
%   in = XPREFERENCES(inputs) sets default values for the input data.
%   Input:  cell array of parameter structures,'inputs'.
%   Output: cell array of parameter structures, with default values set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by xsim
%   Calls qpreferences, phasepreferences, spdpreferences
%   Licensed by Peter D. Drummond, (2023) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sequence = length(inputs);                       % get sequence length
origin = 0;                                      % initial time-origin
for s = 1:sequence                               % loop over sequence 

%  PREFERRED VALUES OF PARAMETERS AND FUNCTIONS

  p = inputs{s};                                 % get parameters for s
  p.version =    xprefer(p,'version',0,'xSPDE4');%% set the version
  p.date    =    xprefer(p,'date',0,datetime);   % set the  date
  p.name =       xprefer(p,'name',0,' ');        % default  name
  p.octave =     xprefer(p,'octave',1,0);        % octave switch
  p.ensembles =  xprefer(p,'ensembles',3,1);     % default ensembles
  p.dimensions = xprefer(p,'dimensions',1,1);    % partial diff. switch
  nd = p.dimensions;
  p.checks =     xprefer(p,'checks',0,[1,zeros(1,nd-1)]);% error checking
  p.tol       =  xprefer(p,'tol',1,1.e-20);      % tolerance in errors
  p.propfactor = xprefer(p,'propfactor',1,@xpropfactor);%%propagator funct.
  if s == 1                                      % if first in sequence?
    fprintf('\n%s: %s\n\n',p.version,p.date);    % print version and date
  else                                           % end if first
    p.checks = inputs{1}.checks;                 % ensure matched checks
  end                                            % end if first
  if isfield(p,'quantum')                        % if quantum switch?
    p = qpreferences(p);                         % get quantum parameters
  else                                           % or else 
    p.quantum = 0;                               % set p.quantum = 0
  end                                            % end if quantum?
  if isfield(p,'phase')                          % if phase-space?
    p = phasepreferences(p);                     % phase-space parameters
  end                                            % end if phase-space
  p.fields =     xprefer(p,'fields',0,{1});      % number of fields
  p.auxfields =  xprefer(p,'auxfields',0,[]);    % auxiliary fields
  p.fieldcells = length(p.fields);               % number of field cells
  p.auxcells   = length(p.auxfields);            % number of auxil. cells
  p.totcells   = p.fieldcells+p.auxcells;        % number of total cells
  p.points     = xprefer(p,'points',1,{[51,35*ones(1,nd-1)]});  
  p.iterations = xprefer(p,'iterations',1,4);    % number of iterations
  p.breed  =     xprefer(p,'breed',1,@Breed);    % default breeding
  p.path   =     xprefer(p,'path',1,@xpath);     % default SDE path
  p.deriv      = xprefer(p,'deriv',0,{@(a,~,~) 0*a});
  p.unoises    = xprefer(p,'unoises',0,[]);      % number of uniform noises
  if isequal(p.unoises ,[])                      % check if uniforms null
    p.noises    = xprefer(p,'noises',p.fieldcells,p.fields);
  else                                           % else uniforms specified
    p.noises    = xprefer(p,'noises',0,[]);      % number of gauss noises
  end                                            % end check uniforms null
  p.fields     = [p.fields,p.auxfields];         % field + auxil. vector
  p.inrandoms  = xprefer(p,'inrandoms',p.fieldcells,p.noises);
  p.ranges     = xprefer(p,'ranges',1,10*ones(1,nd));
  p.steps =      xprefer(p,'steps',p.dimensions,ones(1,nd));
  p.origins =    xprefer(p,'origins',1,[origin,-p.ranges(2:nd)/2]);   
  p.jump   =     xprefer(p,'jump',p.fieldcells,{0});% default jump = 0
  p.rawdata =    xprefer(p,'rawdata',0,0);       % switch for raw data
  p.linear   =   xprefer(p,'linear',p.fieldcells,{0});% linear prop.
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
  p.grid =       xprefer(p,'grid',0,@xgrid);     % handle for setting grid
  p.noisecells  = length(p.noises);
  p.unoisecells = length(p.unoises);
  p.randomcells = length(p.inrandoms);
  p.maxcells = max([p.totcells,p.noisecells+p.unoisecells,p.randomcells]);
  p.methodcells = max([p.fieldcells,p.noisecells+p.unoisecells]);                                          
  if isequal(prod(p.ensembles),1)  
      p.method = xprefer(p,'method',1,@RK4);     % preferred deterministic
    else 
      p.method = xprefer(p,'method',1,@MP);      % preferred stochastic
  end
  p.prop     = xprefer(p,'prop',1,@xprop);       % preferred propagator
                                        
%  CALCULATE CONSTANTS THAT DEPEND ON THE CELL INDEX

  p.nfields      = 0;
  for c = 1:p.fieldcells                        %  index over field cells
    p.nfields = max(length(p.fields{c}),p.nfields);
  end
  p.noisegen =   xprefer(p,'noisegen',1,@xnoise);
  p.randomgen =  xprefer(p,'randomgen',1,@xrandom);
  p.nfilter =    xprefer(p,'nfilter',p.fieldcells,{@xnfilter});
  p.rfilter =    xprefer(p,'rfilter',p.fieldcells,{@xrfilter});
  dx = cell(1,p.maxcells);   a = dx; space = dx;
  ens = p.ensembles(1);
  zs = cell(1,p.noisecells+p.unoisecells);
  for c = 1:p.maxcells                          %  index over all cells
    if length(p.points) < c
        p.points{c} = p.points{c-1};
    end
    p.points{c}(1)  = p.points{1}(1);
    p.points{c}     = [p.points{c},ones(1,nd-length(p.points{c}))];
    dx{c} = p.ranges./max(1,p.points{c}-1);     %  n-th step-size in x
    p.dv{c}  = prod(dx{c}(2:nd));               %  lattice cell volume 
    p.noisefactor{c}  = 1/p.dv{c};
    p.knoisefactor{c} = prod(p.points{c}(2:nd))/p.dv{c};
    space{c}    = p.points{c}(2:p.dimensions);
    p.propagator{c}=1;
    if c <= p.randomcells
      p.inrandoms{c}  = [p.inrandoms{c},0];
      p.d.krandoms{c} = [p.inrandoms{c}(2),space{c},ens];
      p.d.randoms{c}  = [p.inrandoms{c}(1),space{c},ens];
    end
    if c <= p.totcells
      p.d.ca{c}   =  [p.fields{c},space{c},ens];
      p.d.a1{c}   =  [p.fields{c},1,space{c},ens];
      p.d.ft{c}   =  [p.fields{c},p.points{c},ens]*p.spectrum;
      p.d.raw{c}  =  [p.fields{c},p.points{c},ens]*p.rawdata;
      p.d.rc{c}   =  [prod(p.fields{c}),p.points{c},ens];
      a{c}        =  zeros(p.d.ca{c});          %  Initial data for observe
    end    
    if c <= p.noisecells
      p.noises{c}     =  [p.noises{c},0];
      p.d.noises{c}   =  [p.noises{c}(1),space{c},ens];
      p.d.knoises{c}  =  [p.noises{c}(2),space{c},ens];
      zs{c} = zeros([p.noises{c}(1)+p.noises{c}(2),space{c},ens]); 
    end
    if c > p.noisecells && c<= p.unoisecells+p.noisecells
      p.d.noises{c}  =  [p.unoises{c-p.noisecells},space{c},ens];
      zs{c} = zeros(p.d.noises{c}); 
    end
  end
  p.d.r = [ones(1,p.nfields),p.points{1}(2:p.dimensions)]; 
  p.indext = 0;
  p.dx = dx{1};
  p.dk = 2.0*pi./(p.points{1}.*p.dx);           %  n-th step-size in k
  p.kfact =     p.dx/sqrt(2*pi);                %  fft (k) normalization
  p.dimensions = nd;
  p.nspace =     prod(p.points{1}(2:end));      %  count space points
  if p.dimensions > 1
    p = spdpreferences(p);                      %  Partial differential
  else
    p.setboundaries = 0;
  end
  p.dt      = p.ranges(1)/(p.points{1}(1)-1);
  p.nscale  = p.steps(1)*(1+p.checks(1))/p.dt;        %  noise scale factor
  if isequal (p.linear,{[]}) && p.dimensions == 1
      p.prop    =   @(a,~,~) a;
  end
  nc = 1;
  p.checkd(1) = 0;
  p.spacechecks = sum(p.checks(2:nd));
  p.nchk = 1+sum(p.checks);
  for n = 1:p.dimensions
    if p.checks(n)
      nc = nc+1;
      p.checkd(nc) = n;
      p.dxr(n) = p.dx(n)/(2*p.steps(n)); 
    else
      p.dxr(n) = p.dx(n)/(p.steps(n));
    end
  end

%  CALCULATE GENERAL CONSTANTS, 

  p.d.a       =  p.d.ca{1};                     %  Dimensions of first cell
  origin      =  p.origins(1) + p.ranges(1);    %  Origin for next sequence
  p.errors    =  3;                             %  Maximum error field
  p.dw        =  2*pi/(p.ranges(1)+p.dt);
  p.fsc       =  sqrt(2*pi)/p.dw;               %  Fourier scale factor
  p0          =  (0:p.points{1}(1)-1);
  p.plotpts   =  (p0-floor((p.points{1}(1))/2));
  p.w         =  p.plotpts*p.dw;  
  nspace      =  prod(space{1});
  time        =  p.points{1}(1);
  p.xc{1}     =  p.origins(1) + p.dx(1)*p0;     %   1st x-coords
  p.kc{1}     =  p.w;                           %   1st k-coords
  if p.dimensions > 1
      p = p.grid(0,p);                          %  get grid-points
  end
  p.dtr = 0; 
  p.t = p.origins(1);
  [~,s_ord,det_ord,ip,dovect,docell] = p.method(a(1:p.fieldcells),zs,p);
  p.order = xprefer(p,'order',1,0);
  if p.order == -1
    if prod(p.ensembles)>1 
      p.order = s_ord;
    else
      p.order = det_ord;
    end
  end
  p.ipsteps = xprefer(p,'ipsteps',1,ip);
  if ~isequal(p.methodcells,1) && isequal(docell,0)
      error('Cell method needed, total cells = %d\n',p.methodcells);
  end
  if ~isequal(p.fields,1) && isequal(dovect,0)
      error('Vector method needed, total fields = %d\n',p.fields);
  end
 
%  GET PARAMETERS FOR THE OUTPUT AVERAGES

  if ~isfield(p,'observe')
    p.observe =  {@(varargin) varargin{1}};     %   set default observe
  else
    p.observe =  xmcell(p.observe);             %   make observe a cell
  end
  p.expect    =  xprefer(p,'expect',0,{[]});    %   set default expect
  for  e = 1:length(p.expect)
    if ~isequal(p.expect{e},[]) && p.thresholdw %  make weighted expect
      p.observe{e} = @(a,om,p) p.expect{e}(a,p).*exp(om)./mean(exp(om));                   
    end
  end
  p.averages = length(p.observe);
  p.scatters  =  xprefer(p,'scatters',p.averages,{0});
  p.binranges =  xprefer(p,'binranges',p.averages,{{}});%  bin range
  p.transforms = xprefer(p,'transforms',p.averages,{zeros(1,nd+1)});
  p.mincount  =  xprefer(p,'mincount',1,10);
  p.transforms(end+1) = {0};  
  O = cell(1,p.averages);
  p.d.av = cell(1,p.averages);
  for n = 1:p.averages  
    p.spectrum =   p.spectrum || p.transforms{n}(1) ;
    for j = 1:p.dimensions
      if ~p.transforms{n}(j)
        p.xk{n}{j}  = p.origins(j)+(0:p.points{1}(j)-1)*p.dx(j);
      else
        start = -p.points{1}(j)/2;
        p.xk{n}{j}  = (start:start+p.points{1}(j)-1)*p.dk(j);
     end
    end
  end
  
  fprintf('Sequence %d, %s\n',s,p.name);
  fprintf('Calculating %d average(s)\n',p.averages);
 
%  CHECK OBSERVE DIMENSIONS FROM TEST INPUTS

  for n = 1:p.averages                          %  loop over averages
    O{n} = [];                                  %  set observe to void
    if ~isempty(p.observe{n})                   %  if observe defined
      O{n} = p.observe{n}(a{:},p);              %  return observe output
      sz = size(O{n});                          %  size of observe output
      szl = sz(end);                            %  get last observe size
      if szl > ens && rem(szl,ens) == 0         %  check if last size > ens
          sz(end)=szl/ens;                      %  factor: stop mismatch
          sz(end+1)=ens;                         %#ok<AGROW>
      end                                       %  end check last size
      if nd > 1                                 %  if space dimension
        lines = sz(1);                          %  get number of lines
        if prod(sz(2:end-1)) == nspace
          spaceobs = space{1};
          ens1 = sz(end);
        else
          spaceobs = sz(2:end-1);               %  get space points
          ens1 = sz(end);
        end
      else                                      %  if no space dimension
        lines = sz(1:end-1);                    %  get number of lines
        spaceobs = [];                          %  get space points
        ens1 = sz(end);
      end                                       %  end if space
      p.d.obs{n} = [prod(lines),1,prod(spaceobs),ens1];
      spacebins = spaceobs;

%  CALCULATE SCATTER AND PROBABILITY SIZES

      p.scatters{n} = min(p.scatters{n},sz(end));
      if p.scatters{n} > 1                      %   check for scatters 
        lines = lines*p.scatters{n};            %   modify lines
      end                                       %   end check for scatters
      p.d.bins{n} = [];                         %   Initial dimension
      binranges  = p.binranges{n};              %   n-th bin ranges
      p.bins{n} = length(binranges);            %   n-th bin numbers
      if p.bins{n} > 0                          %   If n-th bins set
        if p.verbose > 0                        %   If verbose printing
          fprintf ('\n%d#%d is a %d D probability',s,n,p.bins{n});
        end                                     %   end if verbose printing
        dbin = zeros(1,p.bins{n});              %   dimension of bins
        abin = 1;                               %   Initial area of bin
        for m = 1:p.bins{n}                     %   Loop over bins 
          d = length(binranges{m})-1;           %   Set dimensions of bin 
          p.oc{n}{m} = (binranges{m}(2:d+1)+binranges{m}(1:d))/2;
          p.do{n}{m} = (binranges{m}(2:d+1)-binranges{m}(1:d));
          dbin(m) = d;                          %   Dimension of n-th bin
          abin = abin*p.do{n}{m}(1);            %   Area of n-th bin
          p.xk{n}{nd+m}  = p.oc{n}{m};          %   Make a new axis vector
        end                                     %   End loop over length 
        p.d.bins{n} = dbin;                     %   Dimensions of bin
        p.a.bins{n} = abin;                     %   Area of bin
        spacebins = [spaceobs,dbin];            %   Dimension of space+bins
        lines = 1;                              %   Lines = 1: probability
      end                                       %   End if probability bins
      p.d.av{n} = [lines,time,spacebins];       %   Dimension of averages
      O{n} = zeros(p.d.av{n});                  %   Observe test array cell
    end                                         %   End if observe present
  end                                           %   End loop over averages
  
%  CALCULATE OUTPUT PARAMETERS AND SIZES
    
  p.output =   xprefer(p,'output',0,{[]});
  for n = 1:p.averages
    if n > length(p.output) || isequal(p.output{n},[])
      p.output{n} = @(o,p) o{n};
    end 
  end 
  p.outputs =     length(p.output);
  p.olabels =    xprefer(p,'olabels',p.outputs,{''});%% set the labels
  axes{1} = num2cell(zeros(1,p.dimensions+10));%  initialize default axes
  p.axes =        xprefer(p,'axes',p.outputs,axes);  
  fprintf('Calculating %d output(s)\n\n',p.outputs);
  p.compare =     xprefer(p,'compare',p.outputs,{[]});
  p.ftransforms = xprefer(p,'ftransforms',p.outputs,p.transforms);
  p.binranges =   xprefer(p,'binranges',p.outputs,{{}});
  p.cutoffs  =   xprefer(p,'cutoffs',p.outputs,{-1});
  for n = 1:p.outputs                           %  get output sizes
    if n>p.averages
      p.xk{n} = p.xk{p.averages};               %  get unset coordinates
      p.scatters{n} = 0;
      p.bins{n} = 0;
    end
    if isequal (p.olabels{n}, [])
        p.olabels{n} = ' ';
    end
    p.spectrum = p.spectrum || p.ftransforms{n}(1);
    p.gpoints{n} = size(p.output{n}(O,p));
    p.d.data{n}  = prod(p.gpoints{n});
    p.gpoints{n}(1+end) = max(3,1+p.nchk);
    p.xk{n}{1}  = p.xc{1};
    for j=1:p.dimensions
      if p.ftransforms{n}(j)
        start = -p.points{1}(j)/2;
        p.xk{n}{j}   = (start:start+p.points{1}(j)-1)*p.dk(j);
      end
    end
  end
  inputs{s} = p;
  if p.verbose > 1
    display(p);
    display(p.d);
  end
end                                             %  end loop on sequence
end                                             %  end xpreferences

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END XPREFERENCES FUNCTION

function kn = xnfilter(w,~)                     %   Default noise filters
    kn = w;                                     %   Default noise filter
end

function kr = xrfilter(w,~)                     %   Default random filters
    kr = w;                                     %   Default input filter 
end
