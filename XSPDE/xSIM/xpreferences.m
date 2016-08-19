function in = xpreferences (in)
%   in = XPREFERENCES(in) sets default values for the xspde input data.
%   Input:  input  structure,'in'.
%   Output: input  structure, with default values set. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
  
%%Unconditional  preference list - can be changed if required

    in.version =    xprefer(in,'version',0,'xSIM2.2');
    in.name =       xprefer(in,'name',0,'');
    in.dimension =  xprefer(in,'dimension',1,1);
    in.fields =     xprefer(in,'fields',2,[1,0]);
    in.defines =    in.fields(2);
    in.fields =     in.fields(1);    
    in.fieldsplus = in.fields+in.defines;
    in.ranges =     xprefer(in,'ranges',in.dimension,10*ones(1,in.dimension));
    in.origin =     xprefer(in,'origin',in.dimension,[0,-in.ranges(2:in.dimension)/2.]);
    in.points =     xprefer(in,'points',in.dimension,[51,35*ones(1,in.dimension-1)]);
    in.noises =     xprefer(in,'noises',0,[in.fields,0]);
    in.randoms =    xprefer(in,'randoms',0,in.noises);
    in.ensembles =  xprefer(in,'ensembles',3,[1,1,1]);
    in.boundaries = xprefer(in,'boundaries',in.dimension,zeros(1,in.dimension));
    in.steps =      xprefer(in,'steps',1,1);
    in.iterations = xprefer(in,'iterations',1,4);
    in.order =      xprefer(in,'order',1,1);
    in.checks  =    xprefer(in,'checks',1,[1,zeros(1,in.dimension-1)]);
    in.errorchecks =xprefer(in,'errorchecks',1,sum(in.checks)+1);       %%number of error-check cycles
    in.ebar  =      xprefer(in,'ebar',1,in.errorchecks>1);
    in.octave =     xprefer(in,'octave',1,exist('OCTAVE_VERSION', 'builtin'));
    in.seed =       xprefer(in,'seed',1,0);
    in.file =       xprefer(in,'file',0,'');
    in.print =      xprefer(in,'print',1,1);
    in.raw   =      xprefer(in,'raw',1,0);
    in.numberaxis = xprefer(in,'numberaxis',1,0);
    in.structD =    xprefer(in,'structD',1,0);
%    in.structD =    xprefer(in,'structD',1,1);    %%Structure derivative
    in.errors =     3;                            %%Number of error fields
    
%%Function  preference list - can be changed if required

    in.initial =    xprefer(in,'initial',0,@xinitial);
    in.transfer =   xprefer(in,'transfer',0,@xtransfer);
    if in.structD 
        in.linear=      xprefer(in,'linear',0,@xStructlinear);%%Deprecated
    else 
        in.linear=      xprefer(in,'linear',0,@xlinear);
    end
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
         
%%Conditional preference list - these preferences depend on other choices

    if ~isfield(in,'observe')                %%does input have no observe?
        in.olabels{1} = 'a';
        in.observe{1} =  @(a,~) real(a);     %%default observe function
        if prod(in.ensembles)>1
            in.olabels{1} = '<a>';
        end 
    else                                     %%else input does have observe
        in.observe = xmakecell(in.observe);
        lc = length (in.observe);
        in.olabels = xcprefer(in,'olabels',lc,{' '});
    end
    in.averages =   xprefer(in,'averages',1,length(in.observe));
    in.transforms = xcprefer(in,'transforms',in.averages,{zeros(1,in.dimension)});
    switch func2str(in.step)                 %%set ipsteps for given method
      case {'xEuler','xRK2'}                 %%single transform methods
       in.ipsteps = xprefer(in,'ipsteps',1,1);
      otherwise                              %%double transform methods
       in.ipsteps = xprefer(in,'ipsteps',1,2);
    end
    
  in.functions=0;                            %%Initial plot number
  if isfield(in,'function')                  %%If plot functions input
      in.functions = length(in.function);    %%Number of plot functions
  end                                        %%End if plot functions 
  if in.functions < in.averages              %%More averages than plots
      in.functions = in.averages;            %%Plots set to averages
      in.function{in.functions} =[];         %%Set last plot function
  end                                        %%End averages vs plot
  in.gpoints = xcprefer(in,'gpoints',in.functions,{[1,in.errors,in.points]});
  for n = 1:in.functions                     %% Loop over graphs
    if  isempty(in.function{n})
      if  n<=in.averages
          in.function{n} = @(o,~) o{n};      %% Return default average
      else
          error ('xSIM error: no function, sequence %d, graph %d\n',s,n);
      end
    end                                      %% End if undefined
  end
    
%%Calculated lattice inputs - parameters needed in both xsim and xgraph

    
    in.dx =  in.ranges./max(1,in.points-1);  %%n-th plotted step in x
    in.dk =  2.0*pi./(in.points.*in.dx);     %%n-th step-size in k
    in.dV  = prod(in.dx(2:in.dimension));    %%lattice cell volume
    in.dK  = prod(in.dk(2:in.dimension));    %%k-space volume
    in.nspace = prod(in.points(2:in.dimension));%%Transverse lattice size
    in.V =   in.dV*in.nspace;                %%lattice volume
    in.K  =  in.dK*in.nspace;                %%k-space volume
    for n = 1:in.dimension                   %%loop over  dimension 
      p = (0:in.points(n)-1);                %% index vector
      in.xc{n} = in.origin(n) + p*in.dx(n);  %%n-th x-axis
      p = mod(p+in.points(n)/2,in.points(n));%%n-th cyclic k-index
      p = p - in.points(n)/2;                %%n-th shifted k-index
      in.kc{n} = p*in.dk(n);                 %%n-th propagation k-coords
    end;                                     %%end loop over dimension
end                                          %%end xpreferences function 


%%Default functions - preferred values

function a0 = xinitial(~,r)                  %% First initial data
    a0 = zeros(r.d.a);                       %% set fields to zero
end

function a0 = xtransfer(~,~,a,~)             %% Subsequent data
    a0 = a;                                  %% set to previous output
end

function da = xda(~,~,r)                     %% Default derivative
    da = zeros(r.d.a);                       %% set fields to zero
end

function define = xdefine(~,~,r)             %% Default define
    define = zeros([r.defines,r.nlattice]);  %% set fields to zero
end

function L = xlinear(r)                      %% Default linear filter
    L = zeros(r.d.a);                        %% Default is zero response
end

function Kn = xnfilter(r)                    %% Default stochastic filters
    Kn = ones(r.noises(2),r.nlattice);       %% Default noise filter
end

function Kr = xrfilter(r)                    %% Default stochastic filters
    Kr = ones(r.randoms(2),r.nlattice);      %% Default input filter
end