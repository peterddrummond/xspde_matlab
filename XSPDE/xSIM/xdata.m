function oc = xdata(a,n,r)  
%   o = XDATA(a,n,r) stores data averages from an equation trajectory.
%   Input is the 'a' field and stored data array, 'oc', at current time.
%   Input index n gives the number of the averaging function.
%   Input parameters in the 'r' structure include observe function handles.
%   Returned array 'oc' is the updated observable array, at current time.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
 
tr = r.transforms{n}(2:r.dimension);           %%space transform switch
if sum(tr)>0                                   %%if transform needed
    a = xgraphicsfft(a,tr,r);                  %%Fourier transform
end                                            %%no transform needed   
o1 = r.observe{n}(a,r);                        %%Get x observable    
if r.ensembles(1) > 0
    o1 = reshape(o1,r.d.obs{n}); 
    o1 = mean(o1,2);                           %%Take the mean
end
oc(:,1,:) = o1(:,:);                         %%Return observable
end                                            %%end function
