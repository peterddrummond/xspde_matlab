function input = xpreferences (input)
%   input = XPREFERENCES(input) chooses default values for xspde data.
%   Input:  cell array of sequential input data structures,'input'
%   Output: cell array of data structures, with default values set.
%   Each simulation has its own data structure: input = {in1,in2,...}
%   Probabilities are multidimensional, and binned versus function output
%   See xSPDE manual for complete documentation
%   Called by: xsim, xgraph
%   Needs:     xprefer, xcprefer
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET UNCONDITIONAL PREFERENCES
%
sequence = length(input);                        %%get sequence length
octave = exist('OCTAVE_VERSION','builtin');      %%finds if octave version
t_0  = 0;                                        %%initial default origin 
for s = 1:sequence                               %%loop over sequence 
    in = input{s};                               %%get input structure
    in.version =    xprefer(in,'version',0,'v3.44',0,0);
    in.name =       xprefer(in,'name',0,'',0,0);
    if isfield(in,'dimensions')                 %% If dimensions present
        in.dimension =  in.dimensions;           %% Set dimension 
    end                                          %% Ensures compatibility
    in.dimension =  xprefer(in,'dimension',1,1,1,10);
    nd      =       in.dimension;
    if isfield(in,'fields')  && length(in.fields) > 1 
        if length(in.fields) == 3 
            in.auxfields = in.fields(3);
        end
        in.backfields = in.fields(2);
        in.fields = in.fields(1)+in.fields(2);
    end                                          %% Ensures compatibility
    in.fields =     xprefer(in,'fields',1,1,1,1000000);
    in.backfields = xprefer(in,'backfields',1,0,0,1000000);
    in.auxfields =  xprefer(in,'auxfields',1,0,0,1000000);
    in.breedw    =  xprefer(in,'breedw',1,0,0,0);
    in.mincount  =  xprefer(in,'mincount',1,10,0,0);
    in.fieldsf =    in.fields-in.backfields;
    in.defines  =   in.auxfields; 
    in.fieldsplus = in.fields+in.defines;
    in.ranges =     xprefer(in,'ranges',nd,10*ones(1,nd),0,Inf);
    in.points =     xprefer(in,'points',nd,[51,35*ones(1,nd-1)],2,Inf);
    in.noises =     xprefer(in,'noises',2,[in.fields,0],0,Inf);
    if isfield(in,'inrandoms')                %% If inrandoms present
        in.randoms =  in.inrandoms;           %% Set randoms to inrandoms
    end                                       %% Ensures compatibility
    if isfield(in,'deriv')                   %% If deriv present
        in.da =  in.deriv;                   %% Set da to deriv
    end                                      %% Ensures compatibility
    in.randoms =    xprefer(in,'randoms',0,in.noises,0,Inf);
    in.ensembles =  xprefer(in,'ensembles',3,[1,1,1],1,Inf);
    in.steps =      xprefer(in,'steps',nd,ones(1,nd),1,Inf);
    in.thresholdw = xprefer(in,'thresholdw',1,0,0,0);
    in.iterations = xprefer(in,'iterations',1,4,1,Inf);
    in.iterfb     = xprefer(in,'iterfb',1,100,1,Inf);
    in.iterproj   = xprefer(in,'iterproj',1,4,1,Inf);
    in.checks  =    xprefer(in,'checks',nd,1,0,1); 
    in.errorchecks =xprefer(in,'errorchecks',1,sum(in.checks)+1,1,2);
    in.antialias =  xprefer(in,'antialias',1,0,0,1);
    in.octave =     xprefer(in,'octave',1,octave,0,5);
    in.seed =       xprefer(in,'seed',1,0,0,0);
    in.file =       xprefer(in,'file',0,'',0,0);
    if isfield(in,'verbose')                     %% If verbose present
        in.print =  in.verbose;                  %% Set print to verbose
    end                                          %% Ensures compatibility
    in.print =      xprefer(in,'print',1,0,0,0);
    if isfield(in,'rawdata')                     %% If verbose present
        in.raw =  in.rawdata;                    %% Set print to verbose
    end                                          %% Ensures compatibility
    in.raw   =      xprefer(in,'raw',1,0,0,0);
    in.numberaxis = xprefer(in,'numberaxis',1,0,0,0);
    in.errors =     3;                           %%Set data+error fields
    if isfield(in,'origins')                     %% If origins present
        in.origin =  in.origins;                 %% Set origin to origins
    end                                          %% Ensures compatibility
    in.origin =     xprefer(in,'origin',nd,[t_0,-in.ranges(2:nd)/2],0,0);
    in.origins =    in.origin;
    in.adapt     =  xprefer(in,'adapt',1,1,0,0);
    in.indext =     0;
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET FUNCTION PREFERENCES
%
    in.initial =    xprefer(in,'initial',0,@xinitial,0,0);
    in.initialfb =  xprefer(in,'initialfb',0,@xinitialfb,0,0);
    in.transfer =   xprefer(in,'transfer',0,@xtransfer,0,0);
    in.linear=      xprefer(in,'linear',0,@xlinear,0,0);
    in.noisegen =   xprefer(in,'noisegen',0,@xgaussnoise,0,0);
    in.randomgen =  xprefer(in,'randomgen',0,@xgaussrandom,0,0);
    in.nfilter =    xprefer(in,'nfilter',0,@xnfilter,0,0);
    in.rfilter =    xprefer(in,'rfilter',0,@xrfilter,0,0);
    in.firstfb =    xprefer(in,'firstfb',0,@xfirstfb,0,0);
    if isfield(in,'method')                      %% If method present
        in.step =  in.method;                    %% Set step to method
    end                                          %% Ensures compatibility
    if in.fieldsf == in.fields
     if prod(in.ensembles) == 1
        in.step  = xprefer(in,'step',0,@xRK4,0,0);
      else
        in.step  = xprefer(in,'step',0,@xMP,0,0);
        in.order = xprefer(in,'order',1,1,0,0);
      end
    else
       in.step  = xprefer(in,'step',0,@xMPfb,0,0);
       in.order = xprefer(in,'order',1,0,0,0);
    end
   
    in.grid =       xprefer(in,'grid',0,@xgrid,0,0);
    in.prop =       xprefer(in,'prop',0,@xprop,0,0);
    in.da   =       xprefer(in,'da',0,@xda,0,0);
    in.define   =   xprefer(in,'define',0,@xdefine,0,0);
    in.propfactor = xprefer(in,'propfactor',0,@xpropfactor,0,0);
    in.boundfun =   xprefer(in,'boundfun',0,@xboundfun,0,0);
    in.boundin  =   xprefer(in,'boundin',0,@xboundin,0,0);
    in.points(2:nd) = 1+(in.points(2:nd)-1).*in.steps(2:nd);
    in.dx =         in.ranges./max(1,in.points-1); %%n-th step-size in x
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET CONDITIONAL PREFERENCES        
%
    switch func2str(in.step)                 %%set ipsteps for given method
    case {'xEuler','xImplicit','Euler','Implicit'}%%transform, order 1
       in.ipsteps = 1;
       in.order = xprefer(in,'order',1,1,0,0);
    case {'xRK2','RK2'}                     %%single transform, order 2
       in.ipsteps = 1;
       in.order = xprefer(in,'order',1,2,0,0);
    case {'xMP','xMPadapt','xMPdefine','MP','MPadapt','MPdefine'} 
       in.ipsteps = 2;
       in.order = xprefer(in,'order',1,2,0,0);
    case {'xRK4','RK4'}                            %%double transform, order 4  
       in.ipsteps = 2;
       in.order   = xprefer(in,'order',1,4,0,0);
        otherwise                            %%unknown methods
       in.ipsteps = xprefer(in,'ipsteps',1,0,0,0);
       in.order   = xprefer(in,'order',1,1,0,0);
    end
    if ~isfield(in,'observe')                %%if input has no observe?
        in.observe{1} =  [];                 %%set empty observe function
    end                                      %%end if input has no observe
    if ~isfield(in,'averages')               %%if input has no averages?
        in.averages = max(1,length(in.observe));%%set default
    end                                      %%end if input has no averages
    in.averages = min(in.averages,max(1,length(in.observe)));
    in.observe = xmakecell(in.observe);      %%make observe into a cell
    in.olabels = xmakecell(in.olabels);      %%make olabels into a cell
    if ~isfield(in,'olabels')                %%if input has no olabels?
        in.olabels{in.averages} = [];        %%set empty olabels
    end                                      %%end if input has no olabels
    for n = 1:in.averages                    %%loop over averages
      if  isempty (in.observe{n})            %%if empty observe       
         in.observe{n} =  @(a,~) real(a);    %%set default observe function
         if isempty (in.olabels{n})          %%if no label  
           in.olabels{n} = 'a';              %%set default olabel
           if prod(in.ensembles)>1           %%if input has ensembles?
             in.olabels{n} = '<a>';          %%set default olabel for means
           end                               %%end check ensembles
         end                                 %%end if no label 
      end                                    %%end if empty observe
    end                                      %%end loop over averages
    in.boundaries = xcprefer(in,'boundaries',nd,{zeros(in.fields,2)});
    in.boundaries{1} = -ones(in.fields,2);   %%Neumann time boundary
    in.transforms = xcprefer(in,'transforms',in.averages,{zeros(1,nd+1)});
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE PLOT FUNCTIONS
%
    in.functions=0;                          %%Initial function number
    if isfield(in,'function')                %%If plot functions input
      in.functions = length(in.function);    %%Number of plot functions
    end                                      %%End if plot functions 
    if in.functions < in.averages            %%If too few functions
      in.functions = in.averages;            %%set functions to averages
      in.function{in.functions} =[];         %%Make cell of plot functions
    end                                      %%End if averages > functions
    axes{1} = num2cell(zeros(1,nd+10));      %%initialize default axes
    in.compare = xcprefer(in,'compare',in.functions,{''});
    in.olabels = xcprefer(in,'olabels',in.functions,{' '});
    in.cutoffs = xcprefer(in,'cutoffs',in.functions,{0});
    in.axes =    xcprefer(in,'axes',in.functions,axes);
    in.ftransforms = xcprefer(in,'transforms',in.functions,in.transforms);
    in.scatters =   xcprefer(in,'scatters',in.functions,{0});
    in.binranges = xcprefer(in,'binranges',in.functions,{{}});
    for n = 1:in.functions                   %% Loop over functions
      if  isempty(in.function{n})            %% If no function
        if  n <= in.averages                 %% If average exists
          in.function{n} = @(o,~) o{n};      %% Return default average
        else                                 %% Else return error
          error ('xSIM error: no function, sequence %d, graph %d\n',s,n);
        end                                  %% End if n <= averages
      end                                    %% End if no function
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE SCATTERS, BINS
% 
      binranges  = in.binranges{n};          %% n-th bin ranges matrix
      in.bins{n} = length(binranges);        %% n-th bin numbers
      in.scatters{n} = min(in.scatters{n},in.ensembles(1));
      if in.scatters{n} > 0                  %% If scatters present
        if in.print > 0                      %% If printing
          fprintf ('%d#%d has %d scatters\n',s,n,in.scatters{n});
        end
        if in.ensembles(3) > 1               %% If parallel ensembles
           fprintf('Warning: No scatter-plots if ensembles(3) > 1\n');
           in.scatters{n} = 0;
        end
        if in.bins{n} > 0                    %% If probabilities
           fprintf('Warning: no probabilities with a scatter-plot\n');
           in.bins{n} = 0;                   %% End if probabilities
        end
      end                                    %% End if scatters
      in.d.bins{n} = [];                     %% Initial dimension for bins
      if in.bins{n} > 0                      %% If n-th bins set
         if in.print > 0                     %% If printing
           fprintf ('%d#%d is a %d D probability\n',s,n,in.bins{n});
         end                                 %% end if printing
         dbin = zeros(1,in.bins{n});         %% n-th dimensions for bins
         abin = 1;                           %% Initial area of bin
         for m = 1:in.bins{n}                %% Loop over n-th bin length 
           d = length(binranges{m})-1;       %% Set dimensions of n-th bin
           in.oc{n}{m} = (binranges{m}(2:d+1)+binranges{m}(1:d))/2;
           in.do{n}{m} = (binranges{m}(2:d+1)-binranges{m}(1:d));
           dbin(m) = d;                      %% Dimensions of n-th bin
           abin = abin*in.do{n}{m}(1);       %% Area of n-th bin
         end                                 %% End loop over bin length 
         in.d.bins{n} = dbin;                %% Dimensions of bin results
         in.a.bins{n} = abin;                %% Area of bin
       end
    end                                      %% End loop over functions    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE LATTICE DATA
%
    in.dk =  2.0*pi./(in.points.*in.dx);     %%n-th step-size in k
    in.dkp = pi./((in.points - 1).* in.dx);  %%DCT and DST step
    in.kranges = in.dk.*in.points;           %%ranges in k-space
    in.dv  = prod(in.dx(2:nd));              %%lattice cell volume
    in.dkv = prod(in.dk(2:nd));              %%k-space volume
    in.nspace = prod(in.points(2:nd));       %%Transverse lattice size
    in.v =   in.dv*in.nspace;                %%lattice volume
    in.kv =  in.dkv*in.nspace;               %%k-space volume
    in.points2 = floor(in.points/2);
    for n = 1:nd                             %%loop over  dimension
      p  = (0:in.points(n)-1);               %% index vector
      p1 = (1:in.points(n)-1);               %% truncated index vector
      in.xc{n} = in.origin(n) + p*in.dx(n);  %%n-th x-coords
      in.kc{n} = (p-in.points2(n))*in.dk(n); %%n-th k-coords
      b = in.boundaries{n};                  %%store boundary switch
      type  = b(:,1)+2*b(:,2);               %%get integer switch
      for j = 1:in.fields                    %%loop over fields
        switch type(j)                       %%switch on boundary type
        case -3                              %%Robin boundaries
          in.kcp{j}{n} = p*in.dkp(n);
        case -1                              %%Dirichlet-Robin boundaries
          in.kcp{j}{n} = [0,p1-1/2]*in.dkp(n);
        case 0                               %%periodic boundaries
          %in.kcp{j}{n} = 2*(p-1-in.points2(n))*in.dkp(n);
          %in.kcp{j}{n} = fftshift(in.kcp{j}{n});
          %in.kcp{j}{n}(p1) = fftshift(in.kcp{j}{n}(p1));
          in.kcp{j}{n} = ifftshift(in.kc{n});       %%n-th propagation k-coords
          %testp1 = in.kcp{j}{n};
        case 1                               %%Robin-Dirichlet-boundaries
          in.kcp{j}{n} = [p1-1/2,0]*in.dkp(n);
        case 3                               %%Dirichlet boundaries
          in.kcp{j}{n} = p*in.dkp(n);
        otherwise                            %%invalid boundary type
          error('Invalid IP boundary: [%d,%d]',b(j,1),b(j,2));
        end                                  %%end switch
      end                                    %%end loop over fields
    end                                      %%end loop over dimension
    input{s} = in;                           %%store input structure
    t_0 = in.origin(1) + in.ranges(1);       %%store last time
end                                          %%end sequence loop
end                                          %%end xpreferences function 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XPREFERENCES
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DEFAULT FUNCTION DEFINITIONS
%
function a0 = xinitial(~,r)                  %% First initial data
    a0 = zeros(r.d.a);                       %% set fields to zero
end

function a0 = xinitialfb(~,~,~,r)            %% First initial data
    a0 = zeros(r.d.a);                       %% set fields to zero
end

function [a0,r] = xtransfer(~,r,a,~)         %% Subsequent data
    a0 = a;                                  %% set to previous output
end

function da = xda(~,~,r)                     %% Default derivative
    da = zeros(r.d.a);                       %% set fields to zero
end

function define = xdefine(~,~,r)             %% Default define
    define = zeros([r.d.d]);                 %% set fields to zero
end

function l = xlinear(r)                      %% Default linear filter
    l = zeros(r.d.r);                        %% Default is zero response
end

function kn = xnfilter(w,~)                  %% Default noise filters
    kn = w;                                  %% Default noise filter
end

function kr = xrfilter(w,~)                  %% Default random filters
    kr = w;                                  %% Default input filter 
end

function ap = xfirstfb(a,nc,r)               %% First fb path estimate
    totp = r.steps*nc*(r.points(1)-1)+1;     %% Compute index maximum
    r.d.stored = r.d.fieldsplus;             %% get stored dimensions
    r.d.stored(2) = totp;                    %% update time dimensions
    ap = zeros(r.d.stored);                  %% initialize
    for step = 1:totp                        %% loop over time
        ap(:,step,:) = a;                    %% set to "initial" value
    end                                      %% end loop over time
end                                          %% end function

function b = xboundfun(a,dir,r)              %% Default boundary values
%   b = XBOUNDFUN(a,dir,r) calculates default boundary values
%   Boundaries are set in dimension dir > 0 at upper and lower boundaries
sz = size(a);
dir1 = dir+r.indext;
sz(dir1) = 2;
%% Differentiation direction
if r.t < r.origin(1) ||  r.indext            %% if initial boundary values
    b = zeros(sz);                           %% Set to zeros as default
else                                         %% Not initial, use previous
    b = r.boundval{dir};                     %% Previous boundary values
end 
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END DEFAULT  FUNCTIONS
