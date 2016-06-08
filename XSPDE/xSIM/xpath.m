
function [a,av,raw] = xpath(a,nc,r)
%   [a,data,astore] = XPATH(a,nc,r) solves stochastic equation trajectories.
%   Initial condition is the 'a' field, returned value 'a' is the final field.
%   Input parameters in the 'r' structure including equation function handles.
%   Input 'nc' indicates the check index.
%   Output 'a' is the propagated field, 'o' contains averages and errors, 
%   'raw' is the raw trajectories.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

      %%Initialize stored data

raw = 0;
av = cell(1,r.averages);
if r.raw||r.transformw                         %%if store fields                                     
  raw = zeros(r.d.raw);                        %%initialize storage
end                                            %%end if field stored
for  n = 1:r.averages
    av{n} = zeros(r.d.av{n});
end

      %%Loop over all the time points

for np = 1:r.points(1);                        %%loop until time tmax
  if np > 1                                    %%If NOT first point 
      
      %%Loop over steps within a time point

    for step = 1:r.steps*r.errorchecks;        %%loops over step 
      z = r.noisegen(r);                       %%Calculate  noise
      if (nc<r.errorchecks)&&(step>2*floor(step/2)) %%If low res  & odd
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
        a = reshape(a,r.d.fields);
        raw(:,:,np,:) = a(:,:,1,:);
        a = reshape(a,r.d.a);                  %%reshape fields
  end                                          %%end if store data
  for n = 1:r.averages
     if r.transforms{n}(1) == 0                %%if frequency switch off
        av{n}(:,:,np,:) = xdata(a,n,r);          %%store time-domain data
     end
  end
end;                                           %%end time loop
if r.transformw                                %%if frequency domain
    raw_ft = fft(raw,[],3)*r.kfact(1);         %%take Fourier transform
    for np = 1:r.points(1)                     %%loop until wmax
        aw = reshape(raw_ft(:,:,np,:),r.d.a);  %%flatten field data
        r.w = r.kc{1}(np);                     %%set graphics frequencies
        for n = 1:r.averages
            if r.transforms{n}(1) == 1         %%if frequency switch on
              av{n}(:,:,np,:) = xdata(aw,n,r);
            end
        end
    end                                        %%end loop at wmax
end                                            %%exit frequency domain
end                                            %%end trajectory function

function av = xdata(a,n,r)  
%   av = XDATA(a,n,r) stores data averages from an equation trajectory.
%   Input is the 'a' field of stochastic variables, stored in parallel.
%   Input index n gives the number of the observe function.
%   Input parameters in the 'r' structure include observe function handles.
%   Returned array 'av' is the average observable array, at current time.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
trans = r.transforms{n};                       %%transform switch
if sum(trans(2:r.dimension))>0                 %%if transform needed
    a = xgraphicsfft(a,trans,r);               %%Fourier transform
end                                            %%no transform needed   
o1 = r.observe{n}(a,r);                        %%Get stochastic observable
av = mean(reshape(o1,r.d.obs{n}),2);           %%Take average
end                                            %%end function

function a  =  xgraphicsfft(a,trans,r)            
%   a = XGRAPHICSFFT(a,r,tr) selectively transforms spatial lattice fields.
%   Input is the 'a' field, returned field 'a' is transformed.
%   Input parameters in the 'r' structure including fft phase arrays.
%   Input switch 'tr(i)' = 0 for space domain, = 1 for transform domain.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
a =reshape(a, r.d.fields);                      %%reshape to lattice
for nd = 2:r.dimension                          %%loop over space dimension
    if trans(nd) >0                             %%if FFT required
        a = fft(a,[],2+nd)*r.kfact(nd);         %%take Fourier transform
    end                                         %%end if FFT required
end                                             %%end loop over dimension
a =reshape(a, r.d.a);                           %%reshape to flat array
end                                             %%end function
                