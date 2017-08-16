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
  r.kfacti =    r.dx(1)*r.points(1)/sqrt(2*pi);     %%ifft (w) normalize
  %Note insert correction to phase: exp(iw(t_0+dt/2))
  r.dt  =       r.dx(1)/r.steps;                    %%integration step
  r.dtr  =      r.dt;                               %%reduced step
  r.s.dx =      sqrt(1./r.dv);                      %%Transverse normalize
  r.s.dxt =     r.s.dx*sqrt(r.errorchecks/r.dt);    %%noise normalization
  r.s.dkt =     r.s.dxt*sqrt(r.nspace);             %%k-noise normalization
  r.s.dk  =     r.s.dx*sqrt(r.nspace);              %%input k-noise factor
  r.d.space =   r.points(2:r.dimension);            %%total space dimensions
  r.nspace =    prod(r.d.space);                    %%total space pointss
  r.npoints  =  prod(r.points);                     %%space-time points
  r.d.int =     [r.ensembles(1),1,r.d.space];       %%integration dimensions
  r.nlattice =  r.ensembles(1)*r.nspace;            %%sample*space dimensions
  r.noisetot =  r.noises(1)+r.noises(2);            %%Total noise
  r.d.noises =  [r.noisetot,r.nlattice];            %%Noise dimension
  r.d.r =       [1,r.nlattice];                     %%coordinate dimension
  r.d.a =       [r.fields,r.nlattice];              %%flat field dimension
  r.d.d =       [r.defines,r.nlattice];             %%flat field dimension
  r.d.aplus =   [r.fieldsplus,r.nlattice];          %%flat field dimension
  r.d.fields =  [r.fields,r.d.int];                 %%full field dimension
  r.d.kfields = [r.fields,1,1,r.d.space];           %%broadcastable fields
  r.d.fieldsplus =  [r.fieldsplus,r.d.int];         %%total field 
  r.d.raw =     [r.fieldsplus,r.ensembles(1),r.points];%%raw data  dimension
  r.d.ft =      [r.fieldsplus,r.ensembles(1),r.points(1)-1,r.nspace]; 
  r.infilt =    1;                                  %%default input filter
  r.noisefilt = 1;                                  %%default noise filter
  if r.dimension > 1                                %%if transverse dimensions 
    r = r.grid(r);                                  %%get grid-points
  end                                               %% End check dimension
  
 %%%%%Check the user-supplied functions
  r.t  = r.origin(1);
  w = r.randomgen(r);
  a  = r.initial(w,r);
  if length(a) == 1
            a = a*ones(r.d.a);
  end
  xfcheck('initial',0,a,r.d.a);
  if ~isempty(r.transfer)
      [a,r] = r.transfer(w,r,a,r);
      xfcheck('transfer',0,a,r.d.a);
  end
  d = r.define(a,zeros(r.d.noises),r);
  xfcheck('define',0,d,r.d.d); 
  a(r.fields+1:r.fieldsplus,:) = d;
  da  = r.da(a,r.noisegen(r),r);
  xfcheck('da',0,da,r.d.a);
  lin  = r.linear(r);
  xfcheck('linear',0,lin,r.d.a);
  r.t = 0;
  r.w = 0;
  av = cell(1,r.averages);                          %%number of averages
  for n = 1:r.averages                              %%Loop  over averages
      o = r.observe{n}(a,r);                        %%Size of observe data
      s = xfcheck('observe',n,o,[0,r.nlattice]); 
      r.d.obs{n} = [s(1),r.ensembles(1),1,r.nspace];%%Size of observed data
      r.d.av{n}  = [s(1),1,r.points];               %%Size of averaged data
      av{n} = ones(r.d.av{n});
  end
  for n=1:r.functions
      f = r.function{n}(av,r);
      f = f(:,:);                                   %%Convert to matrix
      s = xfcheck('function',n,f,[0,r.npoints]); 
      r.d.sf{n}    = [s(1),1,r.npoints];
      r.gpoints{n} = [s(1),r.errors,r.points];     %%Add error index
      r.d.data{n}  = [s(1),r.errors,r.npoints]; 
  end
end                                                 %% End function