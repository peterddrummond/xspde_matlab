
function [a,av,raw] = xpathw(a,nc,r)
%   [a,data,astore] = XPATHW(a,nc,r) solves stochastic equation trajectories.
%   Trapezoidal average in time give fourier transform midpoint field values.
%   Initial condition is the 'a' field, returned value 'a' is the final field.
%   Input parameters in the 'r' structure including equation function handles.
%   Input 'nc' indicates the check index.
%   Output 'a' is the propagated field, 'o' contains averages and errors, 
%   'raw' is the raw trajectories.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

      %%Initialize stored data
raw = 0;
av = cell(1,r.averages);
a_ft = zeros(r.d.raw);                      
if r.raw                                       %%if store fields                                     
  raw = zeros(r.d.raw);                        %%initialize storage
end                                            %%end if field stored
for  n = 1:r.averages
    av{n} = zeros(r.d.av{n});
end
a(r.fields+1:r.fieldsplus,:) = r.define(a,zeros(r.d.noises),r); 
totsteps = r.steps*r.errorchecks;
tnorm = 1.0/(r.steps*nc);
astore = a;
t_0 =  r.t;
     
    %%Loop over all the time points

for np = 1:r.points(1);                        %%loop until time tmax 
  if np >1
  r.t = t_0 +(np-2)*r.dt; 
  astore = 0.5*tnorm*a;
    for step = 1:totsteps;                     %%loops over step 
      z = r.noisegen(r);                       %%Calculate  noise
      if (nc<r.errorchecks)&&(step>2*floor(step/2)) %%If low res  & odd
          z1 = z;                              %%Stores noise 
      else                                     %%Else do integration
          if (nc < r.errorchecks)              %%If low res check, even
              z = (z + z1)/2.;                 %%Average noise terms 
          end                                  %%End if low res check
          a = r.step(a,z,r);                   %One step integration
          astore = tnorm*a+astore;
          r.t = r.t+r.dtr;                     %%Increment time
      end                                      %%End if low res  & odd
    end                                        %%End for steps
    astore = astore-0.5*tnorm*a;
  end
  
    %%Store data for the time point
 
 if r.raw                                     %%do fields need storing?
      a =  reshape(a,r.d.fieldsplus);
      raw(:,:,np,:) = a(:,:,1,:);
 end
 astore =  reshape(astore,r.d.fieldsplus);
 a_ft(:,:,np,:) = astore(:,:,1,:);
 
   
   %%compute time-domain observables
   
 for n = 1:r.averages
     if r.transforms{n}(1) == 0                %%if frequency switch off
        av{n}(:,:,np,:) = xdata(a,n,r);        %%store time-domain data
     end
 end
end;                                           %%end time loop
a_ft(:,:,1,:) = 0.5*(a_ft(:,:,1,:)+reshape(a,r.d.fieldsplus));

 
   %%compute frequency-domain observables
   
a_ft = ifft(a_ft,[],3)*r.kfacti;               %%inverse Fourier transform
for np = 1:r.points(1)                         %%loop until wmax
    r.w = r.kc{1}(np);                         %%set graphics frequencies
    astore =  a_ft(:,:,np,:)*exp(1i*(t_0-r.dt/2)*r.w); 
    astore = reshape(astore,r.d.aplus);        %%flatten field data   
    for n = 1:r.averages
        if r.transforms{n}(1) == 1             %%if frequency switch on
            av{n}(:,:,np,:) = xdata(astore,n,r);
        end
    end
end                                            %%end loop at wmax
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
 
a =reshape(a, r.d.fieldsplus);                      %%reshape to lattice
for nd = 2:r.dimension                          %%loop over space dimension
    if trans(nd) >0                             %%if FFT required
        a = fft(a,[],2+nd)*r.kfact(nd);         %%take Fourier transform
    end                                         %%end if FFT required
end                                             %%end loop over dimension
a =reshape(a, r.d.aplus);                           %%reshape to flat array
end                                             %%end function
                