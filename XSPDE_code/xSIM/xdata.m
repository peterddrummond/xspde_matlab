function av = xdata(a,n,r)  
%   av = XDATA(a,n,r) stores data averages, scatter plots and probabilities
%   Input is the 'a' field of stochastic variables, stored in parallel.
%   Input index n gives the number of the observe function.
%   If s=r.scatters{n} > 0, stores a  scatter-plot of first s trajectories.
%   If r.problength{n} > 1, stores probability  density in last dimension.
%   Input parameters in the 'r' structure include observe function handles.
%   Returned array 'av' is the average observable array, at current time.
%   First 'a' dimension is the field index, last dimension the ensemble
%   For probability, the last dimension of 'av' averages indexes over bins
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CALL THE OBSERVE FUNCTION

trans = r.transforms{n};                       %%store transform switch
if sum(trans(2:r.dimension))>0                 %%if transform is needed
    a = xnft(a,trans,r);                       %%make a Fourier transform
end                                            %%end if transform needed
o = r.observe{n}(a,r);                         %%Call the observe function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET PROBABILITY IF NEEDED

pl = r.problength{n};                          %%Probability bins needed
if pl  > 1                                     %%If probability needed
    o = reshape(o,[r.nobs{n},1]);              %%Reshape for binning
    o1r = real(o);                             %%Take the real part
    ub(1,1:pl) =  r.bin{n}(2:pl+1);            %%Upper bin boundary
    lb(1,1:pl) =  r.bin{n}(1:pl);              %%Lower bin boundary
    o=((o1r>lb)&(o1r<ub))./(ub-lb);            %%Probability delta function
end                                            %%End if probability needed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET SCATTER PLOT IF NEEDED

o = reshape(o,[r.d.obs{n},pl]);                %%Reshape for averaging
if r.scatters{n} > 0                           %%Check for scatter plot
    av = zeros(r.scatters{n},r.nspace,1);      %%Initialize scatter plot
    for lines = 1:r.scatters{n}                %%Index over the lines
        av(lines,:,1) = o(1,:,lines);          %%Store data in scatter plot
    end                                        %%End index over the lines
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET AVERAGE IF NEEDED

elseif r.ensembles(1) > 1                      %%Else no scatter plot
    lastind = length(r.d.obs{n});              %%Last index for ensembles
    av = mean(o,lastind);                      %%Take average if possible
else                                           %%No average is possible
    av = o;                                    %%Store current observable
end                                            %%End if scatter plot
end                                            %%End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XDATA FUNCTION