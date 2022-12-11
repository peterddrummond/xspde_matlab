function [a,av,raw] = xpathfb(w,nc,r)
%   [a,av,raw] = XPATHFB(a,nc,r) solves FBSDEs & FBSPDEs by iteration.
%   Initial condition is the 'a' field, returned value 'a' the final field.
%   Initialize requires: a0 = old initial field, a1 = old final field.
%   Fields are a=[ax;ay] where ax propagates forward in time, ay backward
%   Field ay has reverse time order, 'nc' is the check index, 'r' is data.
%   Output 'a' is the final field, 'av' has averages in normal time order, 
%   'avp' contains previous averages, 'raw' are the raw trajectories.
%   xSPDE functions are licensed by Peter D. Drummond, (2020) - see License

raw = zeros(r.d.raw);                            %%initialize raw fields
av = cell(1,2*r.averages);                       %%define average cells                                  
for  n = 1:r.averages                            %%loop over averages  
    av{n} = zeros(r.d.av{n});                    %%initialize averages
end                                              %%end loop over averages
a0=zeros(r.d.a);
a1=zeros(r.d.a);
zstore = storefb(r);                             %%Compute stored noise
totp = r.steps*nc*(r.points(1)-1)+1;             %%Compute index maximum
for iter=1:r.iterfb                              %%loop over iterations
  r.iter = iter;                                 %%store iteration ccounter
  a = r.initialfb(a0,a1,w,r);                    %%Initialize local copy
  if iter == 1                                   %%Check if first iteration
      a0=a;                                      %%Store first boundaries
      ap = r.firstfb(a,nc,r);                    %%First trajectory guess
  end                                            %%End check if first 
  npf=1;                                         %%Initialize time index
  astore = reshape(a,r.d.fieldsplus);            %%Reshape a for storage
  ac(:,1,:) = astore(:,1,:);                     %%Store the initial fields
  npnoise=1;                                     %%Initialize noise index
  r.t = r.origin(1);                             %%Initial time
  for np = 1:r.points(1)                         %%loop until time tmax
    if np > 1                                    %%If after first point 
      for step = 1:r.steps*nc                    %%loops over steps
        z = zstore{npnoise};                     %%Access  noise
        npnoise=npnoise+1;                       %%Increment noise index
        if (nc<r.errorchecks)                    %%If low res
          z = 0.5*(z+zstore{npnoise});           %%Average noise
          npnoise=npnoise+1;                     %%Increment noise index 
        end                                      %%End if low res
        npb=totp-npf+1;                          %%Get time-reversed index
        pr = {ap(:,npb,:),ap(:,npb-1,:);};       %%Old iteration fields
        a = r.step(a,pr,z,r);                    %%One step integration
        astore = reshape(a,r.d.fieldsplus);      %%Reshape fields
        npf=npf+1;                               %%Increment forward index
        ac(:,npf,:) = astore(:,1,:);             %%Store current fields
        r.t = r.t+r.dtr;                         %%Increment time
      end                                        %%End for steps
    end                                          %%End if NOT first point
  raw(:,np,:) = astore(:,1,:);                   %%Store raw fields
  end                                            %%End time loop
  a1 = a;                                        %%Store a1 for iterations
  ap = ac;                                       %%Reset trajectory store
  av(r.averages+1:end) = av(1:r.averages);       %%Reset previous average
  [~,av]=xdatafb(raw,av,iter,r);                 %%Compute all averages
  %if delt<r.fbmaxerror                           %%Check if converged
   % fprintf ('Converged after %d iterations\n',iter);%
   % break
  %end
end                                              %%End iteration loop
end                                              %%end trajectory function

function [zstore] = storefb(r)
%   [zstore] = STOREFB(r) stores random noises 
%   Used for solving FBSDEs by iteration.
%   xSPDE functions are licensed by Peter D. Drummond, (2020) - see License 

total=r.steps*r.errorchecks*(r.points(1)-1);
zstore = cell(1,total);
for npf = 1:total                               %%loops over step
     zstore{npf} = r.noisegen(r);               %%Calculate  noise
end
end

function [del,av] = xdatafb(raw,av,~,r)
%   [av] = XDATAFB(raw,r) stores data from a forward-backward simulation 
%   Used for solving FBSDEs & FBPSDEs by iteration.
%   Calls xdata for each time point, and checks for convergence.
%   xSPDE functions are licensed by Peter D. Drummond, (2020) - see License

for np = 1:r.points(1)                           %%loop until time tmax
  r.t = r.origin(1)+(np-1)*r.dt*r.steps;         %%compute time 
  np1=1+r.points(1)-np;                          %%compute reverse index
  a1(1:r.fieldsf,:) = raw(1:r.fieldsf,np,:);     %%current field forward
  a1(1+r.fieldsf:r.fields,:) = raw(1+r.fieldsf:r.fields,np1,:);%%rev. field
  for n = 1:r.averages                           %%loop over averages  
      av{n}(:,np,:) = xdata(a1,n,r);             %%store time-domain data 
      del = sum(sum(abs(av{n}(:,np,:)-av{n+r.averages}(:,np,:)))); %%error
  end                                            %%end loop over averages
end                                              %%end loop until time tmax
end