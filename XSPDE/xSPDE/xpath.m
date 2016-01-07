
function [a,o,raw] = xpath(a,nc,r)
%   [a,data,astore] = XPATH(a,nc,r) solves stochastic equation trajectories.
%   Initial condition is the 'a' field, returned value 'a' is the final field.
%   Input parameters in the 'r' structure including equation function handles.
%   Input 'nc' indicates the check index.
%   Output 'a' is the propagated field, 'o' contains averages and errors, 
%   'raw' is the raw trajectories.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

dt=r.dt/nc;                                    %%reduced step-size
raw =0;
o = zeros(1,r.points(1),r.nspace,r.graphs);    %%initialize storage
if r.raw||r.transformw                         %%if store fields                                     
  raw = zeros(r.fields,r.nlattice,r.points(1));%%initialize storage
end                                            %%end if field stored                                                                 
for np = 1:r.points(1);                        %%loop until time tmax
  if np > 1                                    %%If NOT first point 
    for step = 1:r.steps*r.errorchecks;        %%loops over step 
      z = r.noisegen(r);                       %%Calculate  noise
      if (nc<r.errorchecks)&&(step>2*floor(step/2))%%If low res  & odd
          z1 = z;                              %%Stores noise 
      else                                     %%Else do integration
          if (nc < r.errorchecks)              %%If low res check, even
              z = (z + z1)/2.;                 %%Average noise terms 
          end                                  %%End if low res check
          a = r.step(a,z,dt,r);                %%One step integration
          r.t = r.t+dt;                        %%Increment time
      end                                      %%End if low res  & odd
    end                                        %%End for steps
  end                                          %%End if NOT first point
  if r.raw||r.transformw                       %%do fields need storing?
        raw(:,:,np) = a;                       %%store fields 
  end                                          %%end if store data
  o(1,np,:,:) = xdata(a,o(1,np,:,:),0,r);      %%store time-domain data
end;                                           %%end time loop
if r.transformw                                %%if frequency domain
    astore = fft(raw,[],3)*r.kfact(1);         %%take Fourier transform
    for np = 1:r.points(1)                     %%loop until wmax
      r.w = r.kc{1}(np);                       %%set graphics frequencies
      o(1,np,:,:) = xdata(astore(:,:,np),o(1,np,:,:),1,r); %%store data
    end                                        %%end loop at wmax
end                                            %%exit frequency domain
end                                            %%end trajectory function

function oc = xdata(a,oc,f,r)  
%   o = XDATA(a,o,f,r) stores data averages from an equation trajectory.
%   Input is the 'a' field and stored data array, 'oc', at current time.
%   Input switch 'f' = 0 for time domain, = 1 for frequency domain averages.
%   Input parameters in the 'r' structure include observe function handles.
%   Returned array 'oc' is the updated observable array, at current time.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
tr1 = zeros(1,r.dimension-1);                  %%initial transform switch
for g = 1:r.graphs                             %%loop over graphs
  if r.transforms{g}(1) == f                   %%if frequency switch match
    tr = r.transforms{g}(2:r.dimension);       %%space transform switch
    if sum(tr)>0                               %%if transform needed
      if sum(tr ~= tr1)                        %%if new space transform
         ak = xgraphicsfft(a,r,tr);            %%Fourier transform
         tr1 = tr;                             %%store transform switch
      end                                      %%end if space transform
      o1 = real(r.observe{g}(ak,r));           %%Get k observable
    else                                       %%no transform needed
      o1 = real(r.observe{g}(a,r));            %%Get x observable    
    end                                        %%end if space transform
    if r.ensembles(1) > 0
        o1 = reshape(o1,r.d.int);              %%Reshape
        o1=mean(o1,1);                         %%Take the mean
    end
    oc(1,1,:,g) = reshape(o1,[r.nspace,1]);         %%Return observable
  end                                          %%end if switches match
end                                            %%end graphs loop
end                                            %%end function

function a  =  xgraphicsfft(a,r,tr)            
%   a = XGRAPHICSFFT(a,r,tr) selectively transforms spatial lattice fields.
%   Input is the 'a' field, returned field 'a' is transformed.
%   Input parameters in the 'r' structure including fft phase arrays.
%   Input switch 'tr(i)' = 0 for space domain, = 1 for transform domain.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt
 
a =reshape(a, r.d.ft);                          %%reshape to lattice
dmax = ndims(a);                                %%get a dimension
for nd = 3:dmax                                 %%loop over dimension
    if tr(nd-2) >0                              %%if FFT required
        a = fft(a,[],nd)*r.kfact(nd-1);         %%take Fourier transform
    end                                         %%end if FFT required
end                                             %%end loop over dimension
a =reshape(a, r.d.a);                           %%reshape to flat array
end                                             %%end function
                