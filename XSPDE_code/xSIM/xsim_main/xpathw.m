function [a,av,raw] = xpathw(a,nc,p)
%   [a,data,astore] = XPATHW(a,nc,p) solves PSDE trajectories.
%   Use for cases where a fourier transform is required in time   
%   Trapezoidal averages in time give midpoint fields for spectra.
%   Initial condition is  'a' field.
%   Input parameters in 'r' structure including equation function handles.
%   Input 'nc' indicates the check index.
%   Output 'a' is the propagated field, 'av' contains averages, 
%   'raw' is for the raw trajectories.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License

%%Initialize stored data

raw = 0;
av = cell(1,p.averages);
aft = zeros(p.d.ft);                             %%initialize FFT store                    
if p.raw                                         %%if store raw fields                                     
  raw = zeros(p.d.raw);                          %%initialize storage
end                                              %%end if field stored
for  n = 1:p.averages
    av{n} = zeros(p.d.av{n});                    %initialize  data
end
nt = p.points(1);
np1 = 1;
nt1 = 1+(nt-1)*p.steps;
tnorm = 0.5/(nc);
t_0 = p.origin(1)-p.dx(1)/(2*p.steps);           %%virtual fft origin
ad = p.define(a,zeros(p.d.noises),p);
 
    %%Loop over all the time points

for np = 1:nt                                    %%loop until time tmax
   if np > 1                                     %%If not first step
     for step = 1:p.steps                        %%loop over steps
       ad = zeros(p.d.d);                        %%initialize the defines
       astore = zeros(p.d.aplus);                %%initialize  stored field
       np1 = np1 + 1;                            %%Count the steps
       for checks = 1:p.errorchecks              %%loop over checks
          z = p.noisegen(p);                     %%Calculate  noise
          if (nc<p.errorchecks)&&(checks == 1)   %%If low res & odd
            z1 = z;                              %%Stores noise but no step 
          else                                   %%Else do integration
            if (nc < p.errorchecks)              %%If low res check, even
              z = (z + z1)/2.;                   %%Average noise terms 
            end                                  %%End if low res check
            if p.defines                         %%if defines needed
              ad = p.define(a,z,p);              %%calculate new defines
            end                                  %%end if defines needed
            astore = astore + tnorm*[a;ad];      %%compute stored field
            a  = p.step(a,z,p);                  %%One step integration
            if p.defines                         %%if defines needed
              ad = p.define(a,z,p);              %%calculate new defines
            end                                  %%end if defines needed
             astore = astore + tnorm*[a;ad];     %%compute stored field
            p.t = p.t+p.dtr;                     %%Increment the time            
          end                                    %%End if low res & odd
       end                                       %%End for checks
       aft(:,np1,:) = astore;                    %%field for transforms
     end                                         %%End for steps
     else                                        %%First time point
     aft(:,1,:) = [a;ad];                        %%store the input field
   end                                           %%End if np >1
   
    %%Store data for the time point  
  
   for n = 1:p.averages                          %%averages loop
     if p.transforms{n}(1) == 0                  %%if frequency switch off
        av{n}(:,np,:) = xdata(a,n,p);            %%store time-domain data
     end                                         %%end if frequency switch
   end                                           %%end averages loop
end                                              %%end np time loop
 
   %%compute frequency-domain observables

storef(:,1,:) = [a;ad];                          %%store the end point
aft(:,1,:) = (aft(:,1,:)+storef(:,1,:))/2;       %%store the cyclic average
if p.raw                                         %%Check if raw needed     
   raw = aft;                                    %%Store raw data            
end                                              %%End if raw
aft = ifft(aft,[],2)*p.kfacti;                   %%Fourier transform
aft = fftshift(aft,2);                           %%shift origin to center
for np = 1:nt                                    %%loop over frequency
    p.w = p.kc{1}(np);                           %%set graphics frequency
    npx = np+floor((nt1-nt)/2);                  %%set frequency point
    astore = aft(:,npx,:)*exp(1i*t_0*p.w); 
    astore = reshape(astore,p.d.aplus);          %%reshape field data
    for n = 1:p.averages
      if p.transforms{n}(1) == 1                 %%if frequency switch on
        av{n}(:,np,:) = xdata(astore,n,p);
      end
    end
end                                              %%end loop at wmax
end                                              %%end trajectory function