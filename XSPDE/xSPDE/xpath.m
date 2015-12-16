
function [a,o,raw] = xpath(a,nc,r)
%   [a,data,astore] = XPATH(a,nc,r) solves stochastic equation trajectories.
%   Initial condition is the 'a' field, returned value 'a' is the final field.
%   Input parameters in the 'r' structure including equation function handles.
%   Input 'nc' indicates the check index.
%   Output 'a' is the propagated field, 'o' contains averages and errors, 
%   'raw' is the raw trajectories.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

dt=r.dt/nc;                                    %%reduced step-size
raw =0;
o = zeros(1,r.points(1),r.n.space,r.graphs);   %%initialize storage
if r.raw||r.transformw                         %%if store fields                                     
    raw = zeros(r.d.raw);                      %%initialize storage
end                                            %%end if field stored                                                                 
for np = 1:r.points(1);                        %%loop until time tmax
  if np > 1                                    %%If NOT first point 
    for step = 1:r.steps*r.errorchecks;        %%loops over step 
      xi = r.noisegen(r);                      %%Calculate  noise
      if (nc<r.errorchecks)&&(step>2*floor(step/2))%%If low res  & odd
          xi1 = xi;                            %%Stores noise 
      else                                     %%Else do integration
          if (nc < r.errorchecks)              %%If low res check, even
              xi = (xi + xi1)/2.;              %%Add noise terms 
          end                                  %%End if low res check
          a = r.step(a,xi,dt,r);               %%One step integration
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
    astore = ifft(raw,[],3)*r.kfact(1);        %%take Fourier transform
    astore = fftshift(astore,3);               %%shift for graphics
    for np = 1:r.points(1)                     %%loop until wmax
      r.w = r.gk{1}(np);                       %%get current frequency
      o(1,np,:,:) = xdata(astore(:,:,np),o(1,np,:,:),1,r); %%store data
    end                                        %%end loop at wmax
end                                            %%exit frequency domain
end                                            %%end trajectory function


function o = xdata(a,o,f,r)  
%   o = XDATA(a,o,f,r) stores data averages from an equation trajectory.
%   Input is the 'a' field, current observable array 'o'.
%   Returned array 'o' is the updated observable array.
%   Input switch 'f' = 0 for time domain, = 1 for frequency domain averages.
%   Input parameters in the 'r' structure include observe function handles.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt
 
tr1 = [0,0,0];                                 %%initial transform switch
for g = 1:r.graphs                             %%loop over graphs
  if r.transforms{g}(1) == f                   %%if frequency switch match
    tr = r.transforms{g}(2:4);                 %%space transform switch
    if sum(tr)>0                               %%if transform needed
      if sum(tr ~= tr1)                        %%if new space transform
         ak = xgraphicsfft(a,r,tr);            %%Fourier transform
         tr1 = tr;                             %%store transform switch
       end                                     %%end if space transform
       o_raw = r.observe{g}(ak,r);             %%Get k observable
    else                                       %%no transform needed
       o_raw = r.observe{g}(a,r);              %%Get x observable    
    end                                        %%end if space transform
    o_raw = reshape(o_raw,[r.ensembles(1),r.n.space]); %%Reshape
    if r.ensembles(1) > 1                      %%If multiple samples
           o_raw=mean(o_raw,1);                %%Take the mean
    end                                        %%End if multiple
    o_raw = reshape(o_raw,[r.n.space,1]);      %%Return observable
    o(1,1,:,g) = real(o_raw(:,1));             %%Return observable
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
                