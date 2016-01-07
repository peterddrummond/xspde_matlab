function input = xinpreferences (input)
%   input = XINPREFERENCES(input) sets default values for the input data.
%   Input:  input cell array of structures,'input'.
%   Output: 'input' cell array, with default values set. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

if isa(input,'cell')                               %%check if input is cell data
    sequence = length(input);                      %%get sequence length of input
else                                               %%else if input is not a cell
    sequence = 1;                                  %%set sequence length to 1  
    input = {input};                               %%change input to cell
end;                                               %%end check type of input
for s = 1:sequence                                 %%loop over sequence length
    in = input{s};                                 %%get input, sequence s
    
%%Unconditional  preference list - can be changed if required

    in.version =    xprefer(in,'version','xSIM1.05');
    in.name =       xprefer(in,'name','');
    in.dimension =  xprefer(in,'dimension',1);
    in.fields =     xprefer(in,'fields',1);
    in.ranges =     xprefer(in,'ranges',10*ones(1,in.dimension));
    in.origin =     xprefer(in,'origin',[0,-in.ranges(2:in.dimension)/2.]);
    in.points =     xprefer(in,'points',[51,35*ones(1,in.dimension-1)]);
    in.noises =     xprefer(in,'noises',[in.fields,0]);
    in.randoms =    xprefer(in,'randoms',in.noises);
    in.ensembles =  xprefer(in,'ensembles',[1,1,1]);
    in.steps =      xprefer(in,'steps',1);
    in.iterations = xprefer(in,'iterations',4);
    in.order =      xprefer(in,'order',1);
    in.errorchecks= xprefer(in,'errorchecks',2);
    in.seed =       xprefer(in,'seed',0);
    in.file =       xprefer(in,'file','');
    in.print =      xprefer(in,'print',1);
    in.raw   =      xprefer(in,'raw',0);
    
%%Function  preference list - can be changed if required

    in.initial =    xprefer(in,'initial',@xinitial);
    in.transfer =   xprefer(in,'transfer',@xtransfer);
    in.linear=      xprefer(in,'linear',@xlinear);
    in.noisegen =   xprefer(in,'noisegen',@xnoisegen);
    in.randomgen =  xprefer(in,'randomgen',@xrandomgen);
    in.nfilter =    xprefer(in,'nfilter',@xnfilter);
    in.rfilter =    xprefer(in,'rfilter',@xrfilter);
    in.step =       xprefer(in,'step',@xMP);
    in.grid =       xprefer(in,'grid',@xgrid);
    in.prop =       xprefer(in,'prop',@xprop);
    in.da   =       xprefer(in,'da',@xda);
    in.propfactor = xprefer(in,'propfactor',@xpropfactor);
    in.numberaxis = xprefer(in,'numberaxis',0);
       
%%Conditional preference list - these preferences depend on other choices


    if ~isfield(in,'observe')                %%does input have no observe?
        in.olabels = xcprefer(in,'olabels',1,{'a_1'});
    else                                     %%else input does have observe
        lc = length (in.observe);
        in.olabels = xcprefer(in,'olabels',lc,{' '});
    end

    in.observe =    xcprefer(in,'observe',1,{@(a,~) real(a(1,:))});
    in.graphs =     xprefer(in,'graphs',length(in.observe));
    in.transforms = xcprefer(in,'transforms',in.graphs,{[0,0,0,0,0,0,0]});
    if s > 1                                 %%check if sequence > 1, 
        in.ensembles = input{s-1}.ensembles; %%set the invariant inputs
        in.order = input{s-1}.order;
    end
    
    switch func2str(in.step)                 %%set ipsteps for given method
      case {'xEuler','xRK2'}                 %%single transform methods
       in.ipsteps = xprefer(in,'ipsteps',1);
      otherwise                              %%double transform methods
       in.ipsteps = xprefer(in,'ipsteps',2);
    end
    
%%Calculated lattice inputs - parameters needed in both xsim and xgraph

    in.dx =  in.ranges./(in.points-1);       %%n-th plotted step in x
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
   input{s}=in;                              %%store input structure
end
end

%%Default functions - preferred values

function a0 = xinitial(~,r)                  %% First initial data
    a0 = zeros(r.d.a);                       %% set fields to zero
end

function a0 = xtransfer(~,~,a,~)             %% Subsequent data
    a0 = a;                                  %% set to previous output
end

function da = xda(~,~,r)                    %% Default derivative
    da = zeros(r.d.a);                       %% set fields to zero
end

function L = xlinear(~,r)                    %% Default linear filter
    L = zeros(r.d.a);                        %% Default is zero response
end

function Kn = xnfilter(r)                %% Default stochastic filters
    Kn = ones(r.d.noise2);                   %% Default noise filter
end

function Kr = xrfilter(r)                %% Default stochastic filters
    Kr = ones(r.d.random2);                  %% Default input filter
end

%%Change history:
%%v0.74 
%%adds default setting for Initialise, 
%%uses xprefer function
%%uses the x-prefix on xspde functions
%%adds default for `store' to store trajectories
%%adds filter output to xlinear for k-space noise
%%adds separate observe functions in Fourier space
%%adds ipstep input
%%v0.91 change to observe cell array
%%v0.95 improved name space 