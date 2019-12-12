function input = xpreferences (input,rawinput)
%   in = XPREFERENCES(in) sets default values for the xspde input data.
%   Input:  input  structure,'in'.
%   Output: input  structure, with default values set.
%   dx=ranges./(step.*points) in is the internal step, 
%   The external graphics step is dxgraph=ranges./points 
%   xSPDE functions are licensed by Peter D. Drummond, (2017) - see License 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET ANY INPUT FROM RAW DATA
%  
if ~isempty(rawinput)                            %% Check rawinput exists
    rawsequence = length(rawinput);              %% Get raw sequence length
    for s= 1:rawsequence                         %% Loop over raw sequence 
        in = rawinput{s};                        %% Get raw sequence input
        fname = fieldnames(in);                  %% Get raw sequence labels
        if s<= length(input)                     %% Check new input length
          for j = 1:length(fname)                %% Loop over raw labels
            label = fname{j} ;                   %% Set new label = raw
            if ~isfield(input{s},label)          %% If no new label data 
              input{s}.(label) = in.(label);     %% Set new input to old
            end                                  %% End if no label data
          end                                    %% End for loop
        else                                     %% No new input available
          input{s} = rawinput{s};                %% Set new input to raw
        end                                      %% End s < = length
    end                                          %% End sequence loop
end                                              %% End if rawinput ~empty
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET UNCONDITIONAL PREFERENCES
%
sequence = length(input);                        %%get sequence length
for s = 1:sequence                               %%loop over sequence 
    in = input{s};                               %%get input structure
    in.version =    xprefer(in,'version',0,'xSIM3.1');
    in.name =       xprefer(in,'name',0,'');
    in.dimension =  xprefer(in,'dimension',1,1);
    nd      =       in.dimension;
    in.fields =     xprefer(in,'fields',2,[1,0]);
    in.defines =    in.fields(2);
    in.fields =     in.fields(1); 
    in.fieldsplus = in.fields+in.defines;
    in.ranges =     xprefer(in,'ranges',nd,10*ones(1,nd));
    in.points =     xprefer(in,'points',nd,[51,35*ones(1,nd-1)]);
    in.noises =     xprefer(in,'noises',2,[in.fields,0]);
    in.randoms =    xprefer(in,'randoms',0,in.noises);
    in.ensembles =  xprefer(in,'ensembles',3,[1,1,1]);
    in.steps =      xprefer(in,'steps',nd,ones(1,nd));
    in.iterations = xprefer(in,'iterations',1,4);
    in.order =      xprefer(in,'order',1,0);
    in.checks  =    xprefer(in,'checks',nd,1);   %%check timestep not space
    in.errorchecks =xprefer(in,'errorchecks',1,sum(in.checks)+1);%%cycles
    in.octave =     xprefer(in,'octave',1,exist('OCTAVE_VERSION','builtin'));
    in.seed =       xprefer(in,'seed',1,0);
    in.file =       xprefer(in,'file',0,'');
    in.rawdata =    xprefer(in,'rawdata',0,'');
    in.print =      xprefer(in,'print',1,1);
    in.raw   =      xprefer(in,'raw',1,0);
    in.numberaxis = xprefer(in,'numberaxis',1,0);
    in.errors =     3;                           %%Number of error fields
    in.broadcast =  xprefer(in,'broadcast',1,0); %%Use array broadcasting
    in.origin =     xprefer(in,'origin',nd,[0,-in.ranges(2:nd)/2]);
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET FUNCTION PREFERENCES
%
    in.initial =    xprefer(in,'initial',0,@xinitial);
    in.transfer =   xprefer(in,'transfer',0,@xtransfer);
    in.linear=      xprefer(in,'linear',0,@xlinear);
    in.noisegen =   xprefer(in,'noisegen',0,@xgaussnoise);
    in.randomgen =  xprefer(in,'randomgen',0,@xgaussrandom);
    in.nfilter =    xprefer(in,'nfilter',0,@xnfilter);
    in.rfilter =    xprefer(in,'rfilter',0,@xrfilter);
    in.step =       xprefer(in,'step',0,@xRK4);
    in.grid =       xprefer(in,'grid',0,@xgrid);
    in.prop =       xprefer(in,'prop',0,@xprop);
    in.da   =       xprefer(in,'da',0,@xda);
    in.define   =   xprefer(in,'define',0,@xdefine);
    in.propfactor = xprefer(in,'propfactor',0,@xpropfactor);
    in.boundfun =   xprefer(in,'boundfun',0,@xboundfun);
    in.boundin  =   xprefer(in,'boundin',0,@xboundin);
    in.points(2:nd) = 1+(in.points(2:nd)-1).*in.steps(2:nd);
    in.dx =         in.ranges./max(1,in.points-1); %%n-th step-size in x
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET CONDITIONAL PREFERENCES        
%
    switch func2str(in.step)                 %%set ipsteps for given method
      case {'xEuler','xRK2','xImplicit'}     %%single transform methods
       in.ipsteps = xprefer(in,'ipsteps',1,1);
      otherwise                              %%double transform methods
       in.ipsteps = xprefer(in,'ipsteps',1,2);
    end
    if ~isfield(in,'observe')                %%does input have no observe?
        in.olabels{1} = 'a';
        in.observe{1} =  @(a,~) real(a);     %%default observe function
        if prod(in.ensembles)>1
            in.olabels{1} = '<a>';
        end 
    else                                     %%else input does have observe
        in.observe = xmakecell(in.observe);  %%make observe into a cell 
    end
    in.boundaries = xcprefer(in,'boundaries',nd,{zeros(in.fields,2)});
    in.averages =   xprefer(in,'averages',1,length(in.observe));
    in.transforms = xcprefer(in,'transforms',in.averages,{zeros(1,nd+1)});
    in.scatters =   xcprefer(in,'scatters',in.averages,{0});
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE PLOT FUNCTIONS
%
    in.functions=0;                          %%Initial function number
    if isfield(in,'function')                %%If plot functions input
      in.functions = length(in.function);    %%Number of plot functions
      %fns = in.functions
    end                                      %%End if plot functions 
    if in.functions < in.averages            %%More averages than functions
      in.functions = in.averages;            %%set functions to averages
      in.function{in.functions} =[];         %%Set last plot function
    end                                      %%End averages vs functions
    in.probability =   xcprefer(in,'probability',in.functions,{0});
    in.olabels = xcprefer(in,'olabels',in.functions,{' '});
    in.ftransforms =xcprefer(in,'transforms',in.functions,in.transforms);
    for n = 1:in.functions                   %% Loop over graphs
    in.problength{n} = length(in.probability{n}); %%Number of probabilities
    if  isempty(in.function{n})
      if  n<=in.averages
          in.function{n} = @(o,~) o{n};      %% Return default average
      else
          error ('xSIM error: no function, sequence %d, graph %d\n',s,n);
      end
    end                                      %% End if undefined
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE LATTICE DATA

    in.dk =  2.0*pi./(in.points.*in.dx);     %%n-th step-size in k
    in.dv  = prod(in.dx(2:nd));              %%lattice cell volume
    in.dkv = prod(in.dk(2:nd));              %%k-space volume
    in.nspace = prod(in.points(2:nd));       %%Transverse lattice size
    in.v =   in.dv*in.nspace;                %%lattice volume
    in.kv =  in.dkv*in.nspace;               %%k-space volume
    in.points2 = floor(in.points/2);
    in.dbounds = 0;
    for n = 1:nd                             %%loop over  dimension
      p = (0:in.points(n)-1);                %% index vector
      in.xc{n} = in.origin(n) + p*in.dx(n);  %%n-th x-coords
      in.kc{n} = (p-in.points2(n))*in.dk(n); %%n-th k-coords
      in.kcp{n} = ifftshift(in.kc{n});       %%n-th propagation k-coords
    end                                      %%end loop over dimension
    nt2 = floor((in.points(1)-1)/2);
    in.kc{1} = in.dk(1)*(-nt2:-nt2+in.points(1)-1);
    input{s} = in;                           %%get input structure
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

function [a0,r] = xtransfer(~,r,a,~)         %% Subsequent data
    a0 = a;                                  %% set to previous output
end

function da = xda(~,~,r)                     %% Default derivative
    da = zeros(r.d.a);                       %% set fields to zero
end

function define = xdefine(~,~,r)             %% Default define
    define = zeros([r.defines,r.nlattice]);  %% set fields to zero
end

function l = xlinear(r)                      %% Default linear filter
    l = zeros(r.d.a);                        %% Default is zero response
end

function kn = xnfilter(w,~)                  %% Default noise filters
    kn = w;                                  %% Default noise filter
end

function kr = xrfilter(w,~)                  %% Default random filters
    kr = w;                                  %% Default input filter 
end

function ab = xboundin(r)                    %% Initial boundary values
    ab = cell(1,r.dimension);                %% set fields to zero
    shape = [r.fields,1,2,1];
    for dir = 2:r.dimension                  %%loop over space dimension
        shape(2) =  prod(r.d.int(1:dir-1)); 
        shape(4) =  prod(r.d.int(dir+1:end)); 
        ab{dir} = zeros(shape);              %%Set to zeros
    end
end

function ab = xboundfun(~,r)                 %% Default boundary values
    ab = r.boundinvalue;                     %% set to initial boundary 
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END DEFAULT  FUNCTIONS
