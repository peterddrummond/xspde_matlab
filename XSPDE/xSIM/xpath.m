
function [a,o,raw] = xpath(a,nc,r)
%   [a,data,astore] = XPATH(a,nc,r) solves stochastic equation trajectories.
%   Initial condition is the 'a' field, returned value 'a' is the final field.
%   Input parameters in the 'r' structure including equation function handles.
%   Input 'nc' indicates the check index.
%   Output 'a' is the propagated field, 'o' contains averages and errors, 
%   'raw' is the raw trajectories.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

raw = 0;
o = cell(1,r.averages);
%%o = zeros(1,r.points(1),r.nspace,r.averages);  %%initialize storage
if r.raw||r.transformw                         %%if store fields                                     
  raw = zeros(r.fields,r.nlattice,r.points(1));%%initialize storage
end                                            %%end if field stored

      %%Loop over all the time points

for np = 1:r.points(1);                        %%loop until time tmax
  if np > 1                                    %%If NOT first point 
      
      %%Loop over all the steps within a single time point

    for step = 1:r.steps*r.errorchecks;        %%loops over step 
      z = r.noisegen(r);                       %%Calculate  noise
      if (nc<r.errorchecks)&&(step>2*floor(step/2))    %%If low res  & odd
          z1 = z;                              %%Stores noise 
      else                                     %%Else do integration
          if (nc < r.errorchecks)              %%If low res check, even
              z = (z + z1)/2.;                 %%Average noise terms 
          end                                  %%End if low res check
          a = r.step(a,z,r);                   %%One step integration
          r.t = r.t+r.dtr;                     %%Increment time
      end                                      %%End if low res  & odd
    end                                        %%End for steps
  end                                          %%End if NOT first point
  
      %%Store data for the time point
        
  if r.raw||r.transformw                       %%do fields need storing?
        raw(:,:,np) = a;                       %%store fields 
  end                                          %%end if store data
  for n = 1:r.averages
     if r.transforms{n}(1) == 0                %%if frequency switch off
        o{n}(1,np,:,:) = xdata(a,n,r);         %%store time-domain data
     end
  end
end;                                           %%end time loop
if r.transformw                                %%if frequency domain
    astore = fft(raw,[],3)*r.kfact(1);         %%take Fourier transform
    for np = 1:r.points(1)                     %%loop until wmax
      r.w = r.kc{1}(np);                       %%set graphics frequencies
      for n = 1:r.averages
        if r.transforms{n}(1) == 1             %%if frequency switch on
            o{n}(1,np,:,:) = xdata(astore(:,:,np),n,r); %%store data
        end
      end
    end                                        %%end loop at wmax
end                                            %%exit frequency domain
end                                            %%end trajectory function

function oc = xdata(a,n,r)  
%   o = XDATA(a,n,r) stores data averages from an equation trajectory.
%   Input is the 'a' field and stored data array, 'oc', at current time.
%   Input index n gives the number of the averaging function.
%   Input parameters in the 'r' structure include observe function handles.
%   Returned array 'oc' is the updated observable array, at current time.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
tr = r.transforms{n}(2:r.dimension);           %%space transform switch
if sum(tr)>0                                   %%if transform needed

    ak = xgraphicsfft(a,r,tr);                 %%Fourier transform
    o1 = real(r.observe{n}(ak,r));             %%Get k observable
else                                           %%no transform needed   
    o1 = real(r.observe{n}(a,r));              %%Get x observable    
end                                            %%end if space transform
s1 = r.ofields(n);
if r.ensembles(1) > 0
    o1 = reshape(o1,[s1,r.d.int]); 
    o1=mean(o1,2);                             %%Take the mean
end
oc(1,1,:,:) = reshape(o1,[s1,r.nspace]).';     %%Return observable
end                                            %%end function

function a  =  xgraphicsfft(a,r,tr)            
%   a = XGRAPHICSFFT(a,r,tr) selectively transforms spatial lattice fields.
%   Input is the 'a' field, returned field 'a' is transformed.
%   Input parameters in the 'r' structure including fft phase arrays.
%   Input switch 'tr(i)' = 0 for space domain, = 1 for transform domain.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
a =reshape(a, r.d.ft);                          %%reshape to lattice
dmax = ndims(a);                                %%get a dimension
for nd = 3:dmax                                 %%loop over dimension
    if tr(nd-2) >0                              %%if FFT required
        a = fft(a,[],nd)*r.kfact(nd-1);         %%take Fourier transform
    end                                         %%end if FFT required
end                                             %%end loop over dimension
a =reshape(a, r.d.a);                           %%reshape to flat array
end                                             %%end function
                