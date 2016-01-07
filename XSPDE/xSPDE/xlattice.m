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
  r.ncopies = r.ensembles(2)*r.ensembles(3);        %%parallel ensembles
  r.transformw = 0;                                 %%initial w-transform
  for i=1:r.graphs                                  %%check for w-transform
      r.transformw = max(r.transformw,r.transforms{i}(1));
  end                                               %%end for w-transform
  r.kfact = r.dx/sqrt(2*pi);                        %%fft normalization
  r.dt  = r.dx(1)/r.steps;                          %%integration step
  r.s.dx = sqrt(1./r.dV);                           %%Transverse normalization
  r.s.dxt = r.s.dx*sqrt(r.errorchecks/r.dt);        %%noise normalization
  r.s.dkt = r.s.dxt*sqrt(r.nspace);                 %%k-noise normalization
  r.s.dk  = r.s.dx*sqrt(r.nspace);                  %%input k-noise factor
  if r.dimension >1                                 %%if transverse dimensions
      r.d.space = r.points(2:r.dimension);          %%total space dimensions
  else                                              %%no transverse dimensions
      r.d.space = 1;                                %%default case
  end                                               %%end if transverse
  r.d.int = [r.ensembles(1),r.d.space];             %%integration dimensions
  r.nlattice = r.ensembles(1)*r.nspace;             %%sample*space dimensions
  r.d.r = [1,r.nlattice];                           %%coordinate dimension
  r.d.noise2 = [r.noises(2),r.d.int];               %%knoise dimension list
  r.d.random2 = [r.randoms(2),r.d.int];             %%krandom dimension list
  r.d.a = [r.fields,r.nlattice];                    %%field dimension list
  r.d.ft = [r.fields,r.d.int];                      %%FT field dimensions
  r.d.obs = [r.ensembles(1),r.nspace,r.graphs];     %%observe dimension list
  r.d.data = [3,r.points(1),r.nspace,r.graphs];     %%data dimension list
  r.d.raw = [r.d.a,r.points(1)];                    %%storage  dimension
  r.infilt = 1;                                     %%default input filter
  r.noisefilt = 1;                                  %%default noise filter
  if r.dimension > 1                                %%if transverse dimensions 
    r = r.grid(r);                                  %%get grid-points
    if (r.noises(2)+r.randoms(2)) > 0               %%if k-noise required
      Kr =  r.rfilter(r);                           %%get k random filters
      Kn =  r.nfilter(r);                           %%get k noise filters
      r.infilt = reshape(Kr,[r.randoms(2),r.nlattice]); %%reshape input filter
      r.noisefilt = reshape(Kn,[r.noises(2),r.nlattice]);%%reshape filter
    end                                             %%end if k-noise 
  end                                               %% End check dimension 
  latt{s} = r;                                      %% Input for sequence
end                                                 %% End loop over sequence
end                                                 %% End function