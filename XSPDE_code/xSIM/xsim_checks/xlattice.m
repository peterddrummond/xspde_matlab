function r = xlattice(r)          
%   latt = XLATTICE(r) sets lattice parameters, dimensions, FFT factors.
%   Input:  input structure, 'r'.
%   Output: input with lattices data, 'r'. 
%   First dimension is the field index, last dimension is the ensemble
%   Note, s.dx = 1/sqrt (dV) = stochastic normalisation in x-space
%   s.dk = 1/sqrt (dK) = sqrt (V/(2*pi)^d) =  normalisation in k-space
%   For (n), if gp = graphed points vector, np = total space points, then: 
%     r.d.obs{n} = [lines,np,ensembles(1)]
%     r.d.av{n}  = [lines,gp(1),...gp(r.dimension)]
%   Called by: xsim
%   Needs:     xfcheck
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STORE INTEGRATION CONSTANTS
  r.ncopies = r.ensembles(2)*r.ensembles(3);     %%parallel ensembles
  r.transformw = 0;                              %%initial w-transform
  for i=1:r.averages                             %%check for w-transform
      r.transformw = max(r.transformw,r.transforms{i}(1));
  end                                            %%end for w-transform
  r.kfact =     r.dx/sqrt(2*pi);                 %%fft (k) normalization
  r.kfspace =   prod(r.kfact(2:end));            %%total k-space norm
  r.kfacti =    r.dx(1)*r.points(1)/sqrt(2*pi);  %%ifft (w) normalize
  r.dt  =       r.dx(1)/r.steps(1);              %%integration step
  r.dtr  =      r.dt;                            %%reduced step
  r.d.space =   r.points(2:r.dimension);         %%total space dimensions
  r.nspace =    prod(r.d.space);                 %%total space points
  r.npoints  =  prod(r.points);                  %%space-time points
  r.s.dx =      sqrt(1./r.dv);                   %%Transverse normalize
  r.s.dxt =     r.s.dx*sqrt(r.errorchecks/r.dt); %%noise normalization
  r.s.dk  =     sqrt(1./r.dkv);                  %%k-random factor 
  r.s.dkt =     r.s.dk*sqrt(r.errorchecks/r.dt); %%k-noise normalization
  r.nlattice =  r.ensembles(1)*r.nspace;         %%sample*space dimensions
  r.d.lattice =  [r.d.space,r.ensembles(1)];     %%integration dimensions
  r.d.int =     [1,r.d.space,r.ensembles(1)];    %%integration dimensions
  r.noisetot =  r.noises(1)+r.noises(2);         %%Total noise
  r.d.noises =  [r.noisetot,r.d.lattice];        %%Noise dimension
  r.d.r =       [1,r.d.space];                   %%coordinate dimension
  r.d.a =       [r.fields,r.d.lattice];          %%field dimension
  r.d.d =       [r.defines,r.d.lattice];         %%defines dimension
  r.d.aplus =   [r.fieldsplus,r.d.lattice];      %%fieldplus dimension
  r.d.fields =  [r.fields,r.d.int];              %%full field dimension
  r.d.fieldsplus =  [r.fieldsplus,1,r.nspace*r.ensembles(1)];
  r.d.raw =     [r.fieldsplus,r.points,r.ensembles(1)];
  r.d.ft =      [r.fieldsplus,r.points(1),r.nspace*r.ensembles(1)]; 
  r.infilt =    1;                               %%default input filter
  r.noisefilt = 1;                               %%default noise filter
  r.d.graph =   [1,r.d.space,r.ensembles(1)];    %%data dimensions
  r.setboundaries = 0;                           %%default boundaries
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BOUNDARIES AND SPATIAL GRID

  for i = 2:r.dimension                          %%Loop over dimensions 
    if ~isequal(r.boundaries{i},zeros(r.fields,2))
       r.setboundaries = 1;                      %%Setboundaries flagged
    end                                          %%End if any boundary set
  end                                            %%End loop over dimension 
  if r.dimension > 1                             %%If space dimensions 
    r = r.grid(0,r);                             %%get grid-points
  end                                            %% End check dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK USER FUNCTIONS ARE OK

  r.t  = r.origin(1)-1;                          %%Set time < origin
  w = r.randomgen(r);                            %%Set random noise
  for dir=2:r.dimension
       r.boundval{dir} = r.boundfun(w,dir,r);    %%Store boundaries
  end
  r.t  = r.origin(1);                            %%Set time = origin
  if r.fields > r.fieldsf
    a  = r.initialfb(zeros(r.d.a),zeros(r.d.a),w,r);
  else 
    a  = r.initial(w,r);
  end
  if ~isempty(r.transfer)
      a = r.transfer(w,r,a,r);
      a = xshape(a,0,r);                         %%reshape fields
  end
  d = r.define(a,zeros(r.d.noises),r);
  d = xshape(d,0,r);                             %%reshape fields
  a = [a;d];
  r.t = r.origin(1)-1;
  r.w = 0;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PROBABILITIES AND AVERAGING

  av = cell(1,r.averages);                       %%Cell for the averages
  for n = 1:r.averages                           %%Loop  over average 
    if ~isempty(r.observe{n})                    %%Check observe existence
      pt = [r.points(1),ones(1,r.dimension-1),r.ensembles(1)];
      o = r.observe{n}(a,r);                     %%Get observe data 
      o = xshape(o,0,r);                         %%Shape the observe data
      s = size(o);                               %%Get the observe size
      pt(2:length(s)) = s(2:end);
      pt = [pt(1:end-1),r.d.bins{n},pt(end)]; 
      if ~isempty(r.d.bins{n})                   %%Check probability bins
          s(1) = 1;                              %%One probability allowed
          nb = prod(r.d.bins{n});
      else
          nb = 1;  
      end                                        %%End check bins
      points2 = prod(pt(2:end-1))*nb;
      r.d.obs{n} = [s(1),points2,pt(end)];       %%Observe points{n}
      r.nobs{n} = prod(r.d.obs{n});              %%Total observe size
      if r.scatters{n} > 0                       %%Check if any scatters
              s(1) = r.scatters{n}*s(1);         %%Get number of scatters
      end                                        %%End check if scatter 
      r.d.av{n}  = [s(1),pt(1:end-1)];           %%Size of averaged data
      av{n} = ones(r.d.av{n});                   %%Initialize averages
    end                                          %%End check observes
  end                                            %%End loop  over averages
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE FUNCTIONAL OUTPUT
  
  r.indext = 1;
  for n = 1:r.functions
      pt = ones(1,r.dimension);
      f = r.function{n}(av,r);
      f = xshape(f,1,r);                         %%Shape the function data
      s = size(f); 
      pt(1:length(s)-1) = s(2:end);
      npoints=prod(pt);
      r.d.sf{n}    = [s(1)*npoints,1];
      r.gpoints{n} = [s(1),pt,r.errors];         %%Add error index
      r.d.data{n}  = [s(1)*npoints,r.errors]; 
  end
  r.indext = 0;
end                                              %% End function