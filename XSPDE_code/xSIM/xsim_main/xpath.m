function [a,av,raw] = xpath(a,nc,p)
%   [a,data,astore] = XPATH(a,nc,p) solves stochastic equation trajectories.
%   Initial condition is the 'a' field, returned value 'a' is the final field.
%   Input parameters in the 'r' structure including equation function handles.
%   Input 'nc' indicates the check index.
%   Output 'a' is the propagated field, 'o' contains averages and errors, 
%   'raw' is the cell array of raw trajectories.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

      %%Initialize stored data

raw=0;
av = cell(1,p.averages);
if p.raw                                         %%if raw fields are stored                                    
  raw = zeros(p.d.raw);                          %%initialize raw storage
end                                              %%end if raw field stored
for  n = 1:p.averages
    av{n} = zeros(p.d.av{n});
end
z = zeros(p.d.noises);
ad = p.define(a,z,p);
tnorm = 0.5/(nc);
for np = 1:p.points(1)                           %%loop until time tmax
  if np > 1                                      %%If NOT first point
    if p.thresholdw > 0 && np > 2
        [a,p.breedw] = xbreed(a,p);              %Breeding by weights
    end
      
      %%Loop over steps within a time point
 
    for step = 1:p.steps(1)                      %%loops over step
      ad = 0.*ad;
      for checks = 1:p.errorchecks               %%loop over checks
        z1 = z;                                  %%Calculate  noise
        z = p.noisegen(p);                       %%Calculate  noise
        if (nc>=p.errorchecks)||(checks == 2)    %%If ~( low res  & odd)
          if (nc < p.errorchecks)                %%If low res check, even
              z = (z + z1)/2.;                   %%Average noise terms 
          end                                    %%End if low res check
          if p.defines                           %%if defines needed
              ad = ad+tnorm*p.define(a,z,p);     %%calculate new defines
          end                                    %%end if defines needed
          a = p.step(a,z,p);                     %%One step integration          
          p.t = p.t+p.dtr;                       %%Increment time
          if p.defines                           %%if defines needed
              ad = ad+tnorm*p.define(a,z,p);     %%calculate new defines      
          end                                    %%End if defines needed
        end                                      %%End if low res & odd
      end                                        %%End for checks
    end                                          %%End for steps
  end                                            %%End if NOT first point
  
      %%Store data for the time point

  if p.raw
      astore = reshape([a;ad],p.d.fieldsplus);
      raw(:,np,:) = astore(:,1,:);               %#ok<AGROW>
  end                                            %%end if store data
  for n = 1:p.averages
      av{n}(:,np,:) = xdata([a;ad],n,p);         %%store time-domain data
  end
end                                              %%end time loop
end                                              %%end trajectory function