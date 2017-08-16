
function [a,av,raw] = xpathw(a,nc,r)
%   [a,data,astore] = XPATHW(a,nc,r) solves PSDE trajectories.
%   Use for cases where a fourier transform is required in time   
%   Trapezoidal averages in time give midpoint field values.
%   Initial condition is  'a' field.
%   Input parameters in 'r' structure including equation function handles.
%   Input 'nc' indicates the check index.
%   Output 'a' is the propagated field, 'av' contains averages, 
%   'raw' is the raw trajectories.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

      %%Initialize stored data
raw = 0;
av = cell(1,r.averages);
a_ft = zeros(r.d.ft);                      
if r.raw                                       %%if store fields                                     
  raw = zeros(r.d.raw);                        %%initialize storage
end                                            %%end if field stored
for  n = 1:r.averages
    av{n} = zeros(r.d.av{n});
end
a(r.fields+1:r.fieldsplus,:) = r.define(a,zeros(r.d.noises),r);
totsteps = r.steps*r.errorchecks;
tnorm = 0.5/(r.steps*nc);
nt = r.points(1);
t_0 = r.origin(1);
     
    %%Loop over all the time points

for np = 1:nt;                        %%loop until time tmax 
  if np >1
    astore = zeros(r.d.aplus);                 %%initialize stored field 
    for step = 1:totsteps;                     %%loops over step 
      z = r.noisegen(r);                       %%Calculate  noise
      if (nc<r.errorchecks)&&(step>2*floor(step/2)) %%If low res  & odd
          z1 = z;                              %%Stores noise, no step 
      else                                     %%Else do integration
          if (nc < r.errorchecks)              %%If low res check, even
              z = (z + z1)/2.;                 %%Average noise terms 
          end                                  %%End if low res check
          a=[a(1:r.fields,:);r.define(a,z,r)]; %%get the defined fields  
          astore = tnorm*a+astore;             %%average stored field 
          a = r.step(a,z,r);                   %One step integration
          astore = tnorm*a+astore;             %%average stored field 
          r.t = r.t+r.dtr;                     %%Increment time
      end                                      %%End if low res  & odd
    end                                        %%End for steps
    astore =  reshape(astore,r.d.fieldsplus);  %%reshape stored field 
    a_ft(:,:,np-1,:) = astore(:,:,1,:);        %%store field for Fourier
  end                                          %%End if np >1
  
    %%Store data for the time point
    
 if r.raw         
      astore = reshape(a,r.d.fieldsplus);
      raw(:,:,np,:) = astore(:,:,1,:);         %%Save raw trajectories
 end                                           %%End if raw
 for n = 1:r.averages
     if r.transforms{n}(1) == 0                %%if frequency switch off
        av{n}(:,:,np,:) = xdata(a,n,r);        %%store time-domain data
     end                                       %%end if frequency switch
 end                                           %%end averages loop
end                                            %%end time loop
 
   %%compute frequency-domain observables

a_ft = ifft(a_ft,[],3)*r.kfacti;               %%inverse Fourier transform
a_ft = fftshift(a_ft,3);                       %%shift origin to center
a_ft(:,:,nt,:) = a_ft(:,:,1,:);                %%symmetrize for plotting
for np = 1:nt                                  %%loop over frequency
    r.w = r.kc{1}(np);                         %%set graphics frequency
    astore =  a_ft(:,:,np,:)*exp(1i*(t_0+r.dt/2)*r.w); 
    astore = reshape(astore,r.d.aplus);        %%flatten field data
    for n = 1:r.averages
      if r.transforms{n}(1) == 1               %%if frequency switch on
        av{n}(:,:,np,:) = xdata(astore,n,r);
      end
    end
end                                            %%end loop at wmax
end                                            %%end trajectory function