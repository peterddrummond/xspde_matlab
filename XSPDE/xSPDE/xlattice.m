function latt = xlattice(input)          
%   latt = XLATTICE(input) sets lattice parameters, dimensions, FFT factors.
%   Input:  input cell array, 'input'.
%   Output: lattice cell array, 'latt'. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

sequence = length(input);
latt=cell(sequence);                                %%lattice cell
for s = 1:sequence                                  %%loop over sequence
  in = input{s};                                    %%extract sequence s
  r = in;                                           %%copy to r struct
  r.n.ensemble = r.ensembles(2)*r.ensembles(3);     %%parallel ensembles
  r.transformw = 0;
  for i=1:r.graphs
        r.transformw = max(r.transformw,r.transforms{i}(1));
  end 
  
  kph={1};                                          %%k phase-factors  FFT
  xph={1};                                          %%x phase-factors  FFT
  s2pi = sqrt(2*pi);                                %%square root 2 pi
  if r.transformw
    r.w = r.gk{1};                                   %%set frequency data
    r.wtph=s2pi*exp(-1i*0.5*r.kr(1)*r.xc{1})/r.dk(1);%%wt phase-factors
  end
  for n= 2:r.dimension                              %%loop over  dimension
    minus=ones(1,prod(r.points(2:n-1)));            %%low dimension  points
    plus=ones(1,prod(r.points(n+1:r.dimension)));   %%high dimension  points 
    xph{n}=exp(1i*0.5*r.kr(n)*(r.xc{n}-r.origin(n)));%% x phase-factor
    kph{n}=r.dx(n)*exp(-1i*r.gk{n}*r.origin(n))/s2pi;%% k phase-factor
    xph{n}=kron(plus,kron(xph{n},minus));           %% x phase-factor fill
    kph{n}=kron(plus,kron(kph{n},minus));           %% k phase-factor fill
  end;                                              %%end loop dimension
  r.dt  = r.dx(1)/r.steps;                          %%integration step
  r.s.dx = sqrt(1./r.dV);                           %%Transverse cell
  r.s.dxt = r.s.dx*sqrt(r.errorchecks/r.dt);        %%noise normalization
  r.n.space = prod(r.points(2:r.dimension));        %%Transverse lattice 
  r.s.dkt = r.s.dxt*sqrt(r.n.space);                %%k-noise normalization
  r.s.dk  = r.s.dx*sqrt(r.n.space);                 %%input k-noise factor
  if r.dimension >1                                 %%if transverse dimensions
      r.d.space = r.points(2:r.dimension);          %%total space dimensions
  else                                              %%no transverse dimensions
      r.d.space = 1;                                %%default case
  end                                               %%end if transverse
  r.d.int = [r.ensembles(1),r.d.space];             %%integration dimensions
  r.n.lattice = r.ensembles(1)*r.n.space;           %%sample*space dimensions
  r.d.r = [1,r.n.lattice];                          %%coordinate dimension
  r.d.noise2 = [r.noises(2),r.d.int];               %%knoise dimension list
  r.d.random2 = [r.randoms(2),r.d.int];             %%krandom dimension list
  r.d.a = [r.fields,r.n.lattice];                   %%field dimension list
  r.d.ft = [r.fields,r.d.int];                      %%FT field dimensions
  r.d.obs = [r.ensembles(1),r.n.space,r.graphs];    %%observe dimension list
  r.d.data = [3,r.points(1),r.n.space,r.graphs];    %%data dimension list
  r.d.raw = [r.d.a,r.points(1)];                    %%storage  dimension
  if r.dimension > 1                                %%if transverse dimensions 
    r = r.grid(r);                                  %%get grid-points
    if (r.noises(2)+r.randoms(2)) > 0               %%if k-noise required
      Kr =  r.rfilter(r);                           %%get k random filters
      Kn =  r.nfilter(r);                           %%get k noise filters
      r.infilt = reshape(Kr,[r.randoms(2),r.n.lattice]); %%reshape input filter
      r.noisefilt = reshape(Kn,[r.noises(2),r.n.lattice]); %%reshape noise filter
    end                                             %%end if k-noise 
    fmat = ones(r.fields,r.ensembles(1));           %%Matrix: field, ensembles
    for n = 2:r.dimension                           %%loop over  dimension 
      r.pre{n} = reshape(kron(xph{n},fmat),r.d.ft); %%Graphics FT prefactor
      r.post{n}= reshape(kron(kph{n},fmat),r.d.ft); %%Graphics FT postfactor
    end;                                            %%end loop dimension
  end                                               %% End check dimension 
  latt{s} = r;                                      %%input for sequence
end                                                 %%End loop over sequence
end                                                 %%End function
 
%%Changes in 0.41
%%Increased decimal places in heading, add sample number to heading
%%Changes in 0.51
%%Remove graphics data, dimension data in al struct, full dimensions
%%Changes in 0.6
%%Add sequence processing, changed definition of ftc cell array
%%Changes in 0.61
%%Uses dt,dw in standard normalization of dw^2~dt
%%Changes in 0.62
%%Fix bug with mutiple fields in one dimension
%%Changes in 0.7 
%%shorten the v0.62 code, return L default to preferences
%%Changes in 0.71: improve definition of K lattice
%%Changes in 0.72: reduce number of dimensions in dma.
%%Changes in 0.8: add ipsteps logic, modify k normalization