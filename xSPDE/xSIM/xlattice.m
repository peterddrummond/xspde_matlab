function r = xlattice(r)          
%   latt = XLATTICE(in) sets lattice parameters, dimensions, FFT factors.
%   Input:  input structure, 'r'.
%   Output: input with lattices data, 'r'. 
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STORE INTEGRATION CONSTANTS
  
  r.ncopies = r.ensembles(2)*r.ensembles(3);      %%parallel ensembles
  r.transformw = 0;                               %%initial w-transform
  for i=1:r.averages                              %%check for w-transform
      r.transformw = max(r.transformw,r.transforms{i}(1));
  end                                             %%end for w-transform
  r.kfact =     r.dx/sqrt(2*pi);                  %%fft (k) normalization
  r.kfacti =    r.dx(1)*r.points(1)/sqrt(2*pi);   %%ifft (w) normalize
  %Note insert correction to phase: exp(iw(t_0+dt/2))
  r.dt  =       r.dx(1)/r.steps(1);               %%integration step
  r.dtr  =      r.dt;                             %%reduced step
  r.s.dx =      sqrt(1./r.dv);                    %%Transverse normalize
  r.s.dxt =     r.s.dx*sqrt(r.errorchecks/r.dt);  %%noise normalization
  r.s.dkt =     r.s.dxt*sqrt(r.nspace);           %%k-noise normalization
  r.s.dk  =     r.s.dx*sqrt(r.nspace);            %%input k-noise factor
  r.d.space =   r.points(2:r.dimension);          %%total space dimensions
  r.nspace =    prod(r.d.space);                  %%total space pointss
  r.npoints  =  prod(r.points);                   %%space-time points
  r.d.int =     [1,r.d.space,r.ensembles(1)];     %%integration dimensions
  r.nlattice =  r.ensembles(1)*r.nspace;          %%sample*space dimensions
  r.noisetot =  r.noises(1)+r.noises(2);          %%Total noise
  r.d.noises =  [r.noisetot,r.nlattice];          %%Noise dimension
  r.d.r =       [1,r.nlattice];                   %%coordinate dimension
  r.d.a =       [r.fields,r.nlattice];            %%flat field dimension
  r.d.d =       [r.defines,r.nlattice];           %%defines dimension
  r.d.aplus =   [r.fieldsplus,r.nlattice];        %%fieldplus dimension
  r.d.fields =  [r.fields,r.d.int];               %%full field dimension
  r.d.fieldsplus =  [r.fieldsplus,r.d.int];       %%total field 
  r.d.raw =     [r.fieldsplus,r.points,r.ensembles(1)];%%raw data dimension
  r.d.ft =      [r.fieldsplus,r.points(1)-1,r.nspace,r.ensembles(1)]; 
  r.infilt =    1;                                %%default input filter
  r.noisefilt = 1;                                %%default noise filter
  r.setboundaries = 0;                            %%default boundaries
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BOUNBARIES AND SPATIAL GRID

  for i = 2:r.dimension                           %%Loop over dimensions 
    if ~isequal(r.boundaries{i},zeros(r.fields,2))%%If any boundary is set
       r.setboundaries = 1;                       %%Setboundaries flagged
    end                                           %%End if any boundary set
  end                                             %%End loop over dimension 
  if r.dimension > 1                              %%If space dimensions 
    r = r.grid(1,r);                              %%get grid-points
  end                                             %% End check dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK USER FUNCTIONS ARE OK

  r.t  = r.origin(1);                             %%Set time at origin
  w = r.randomgen(r);                             %%Set random noise
  r.boundinvalue = r.boundin(r);                  %%Set boundaries
  a  = r.initial(w,r);
  if length(a) == 1
            a = a*ones(r.d.a);
  end
  xfcheck('initial',0,a,r.d.a);                   %%Check initialize func
  if ~isempty(r.transfer)
      [a,r] = r.transfer(w,r,a,r);
      xfcheck('transfer',0,a,r.d.a);              %%Check transfer func
  end
  d = r.define(a,zeros(r.d.noises),r);
  xfcheck('define',0,d,r.d.d);                    %%Check define func
  a = [a;d];
  da  = r.da(a,r.noisegen(r),r);
  xfcheck('da',0,da,r.d.a);                       %%Check derivative func
  lin  = r.linear(r);
  xfcheck('linear',0,lin,r.d.a);                  %%Check linear func
  r.t = 0;
  r.w = 0;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PROBABILITIES AND AVERAGING

  av = cell(1,r.averages);                        %%Cell for the averages
  points=av;                                      %%Cell for the points
  for n = 1:r.averages                            %%Loop  over averages
      o = r.observe{n}(a,r);                      %%Get the observe data
      s = xfcheck('observe',n,o,[0,r.nlattice]);  %%Check the observe data
      r.d.obs{n} = [s(1),r.nspace,r.ensembles(1)];%%Size of observed data
      r.nobs{n} = prod(r.d.obs{n});
      r.scatters{n} = min(r.scatters{n},r.ensembles(1));
      if r.scatters{n} > 0   
        s(1) = r.scatters{n};                     %%Get number of scatters
      end
      points{n}=r.points;
      lp = r.problength{n};                       %%Number of probabilities
      if lp >1                                    %%If probability needed
          r.bin{n}=zeros(1,lp+1);
          p=r.probability{n};
          r.bin{n}(2:lp) = (p(2:lp)+p(1:lp-1))/2; %% bin boundary
          r.bin{n}(1)    = p(1);                  %% lower bin boundary
          r.bin{n}(lp+1) = p(lp);                 %% upper bin boundary  
	      points{n}=   [r.points,lp];
      end 
      r.d.av{n}  = [s(1),points{n},1];            %%Size of averaged data   
      av{n} = ones(r.d.av{n});                    %%Intialize averages
  end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE FUNCTIONAL OUTPUT

  for n=1:r.functions
      f = r.function{n}(av,r);
      f = f(:,:);                                  %%Convert to matrix
      if n>length(points)
          points{n}=r.points;
      end
      npoints=prod(points{n});
      s = xfcheck('function',n,f,[0,npoints]); 
      r.d.sf{n}    = [s(1)*npoints,1];
      r.gpoints{n} = [s(1),points{n},r.errors];    %%Add error index
      r.d.data{n}  = [s(1)*npoints,r.errors]; 
  end
end                                                 %% End function