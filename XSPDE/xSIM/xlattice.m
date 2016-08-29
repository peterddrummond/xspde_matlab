function r = xlattice(r)          
%   latt = XLATTICE(in) sets lattice parameters, dimensions, FFT factors.
%   Input:  input structure, 'r'.
%   Output: input with lattices data, 'r'. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

  r.ncopies = r.ensembles(2)*r.ensembles(3);        %%parallel ensembles
  r.transformw = 0;                                 %%initial w-transform
  
  for i=1:r.averages                                %%check for w-transform
      r.transformw = max(r.transformw,r.transforms{i}(1));
  end                                               %%end for w-transform
  r.kfact =     r.dx/sqrt(2*pi);                    %%fft (k) normalization
  r.kfacti =    r.ranges(1)/sqrt(2*pi);             %%ifft (w) normalize
  %Note insert correction to phase: exp(iw(t_0+dt/2))
  r.dt  =       r.dx(1)/r.steps;                    %%integration step
  r.s.dx =      sqrt(1./r.dv);                      %%Transverse normalize
  r.s.dxt =     r.s.dx*sqrt(r.errorchecks/r.dt);    %%noise normalization
  r.s.dkt =     r.s.dxt*sqrt(r.nspace);             %%k-noise normalization
  r.s.dk  =     r.s.dx*sqrt(r.nspace);              %%input k-noise factor
  r.d.space =   r.points(2:r.dimension);            %%total space dimensions
  r.npoints  =  prod(r.points);                     %%total points
  r.d.int =     [r.ensembles(1),1,r.d.space];       %%integration dimensions
  r.nlattice =  r.ensembles(1)*r.nspace;            %%sample*space dimensions
  r.noisetot =  r.noises(1)+r.noises(2);            %%Total noise
  r.d.noises =  [r.noisetot,r.nlattice];            %%Noise dimension
  r.d.r =       [1,r.nlattice];                     %%coordinate dimension
  r.d.a =       [r.fields,r.nlattice];              %%flat field dimension
  r.d.aplus =   [r.fieldsplus,r.nlattice];          %%flat field dimension
  r.d.fields =  [r.fields,r.d.int];                 %%full field dimension
  r.d.fieldsplus =  [r.fieldsplus,r.d.int];          %%total field 
  r.d.raw =     [r.fieldsplus,r.ensembles(1),r.points]; %%raw data  dimension
  r.infilt =   1;                                   %%default input filter
  r.noisefilt =1;                                   %%default noise filter
  if r.dimension > 1                                %%if transverse dimensions 
    r = r.grid(r);                                  %%get grid-points
    if (r.noises(2)+r.randoms(2)) > 0               %%if k-noise required
      Kr =  r.rfilter(r);                           %%get k random filters
      Kn =  r.nfilter(r);                           %%get k noise filters
      r.infilt = reshape(Kr,[r.randoms(2),r.d.int]); %%reshape input filter
      r.noisefilt = reshape(Kn,[r.noises(2),r.d.int]);%%reshape filter
    end                                             %%end if k-noise 
  end                                               %% End check dimension
  a = zeros(r.d.aplus);
  r.t = 0;
  av = cell(1,r.averages);
  for n = 1:r.averages                              %%Loop  over averages
      so = size(r.observe{n}(a,r));                 %%Size of observe data
      if length(so) == 2 
            ns = so(2)/r.ensembles(1);
      else
            ns =  prod(so(3:end));
      end 
      r.d.obs{n} = [so(1),r.ensembles(1),1,ns];     %%Size of observed data
      r.d.av{n}  = [so(1),1,r.points(1),ns];        %%Size of averaged data
      av{n} = zeros([so(1),1,r.points(1)*ns]);
  end
  for n=1:r.functions
      sf = size(r.function{n}(av,r));               %%Size of function data
      if numel(sf) == 3 &&  sf(3) == r.npoints
               sf =[sf(1),r.points];                %%Expand function size
      end
      npoints      = prod(sf(2:end));
      r.d.sf{n}    = [sf(1),1,npoints];
      r.gpoints{n} = [sf(1),r.errors,sf(2:end)];   %%Add error index
      r.d.data{n}  = [sf(1),r.errors,npoints]; 
  end
end                                                 %% End function