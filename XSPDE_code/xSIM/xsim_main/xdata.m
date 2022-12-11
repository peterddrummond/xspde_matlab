function av = xdata(a,n,p)  
%   av = XDATA(a,n,p) stores data averages, scatter plots or probabilities
%   Input is the 'a' field of stochastic variables, stored in parallel:
%   where the 'a' dimensions are:[field index, space*ensemble index].
%   Input index n gives the number of the observe function.
%   Input struct r gives parameters and observe function handles.
%   Observe dimension is [line index, space*ensemble index].
%   If s=p.scatters{n} > 0, stores a  scatter-plot of first s trajectories.
%   If p.bin{n} is not empty, stores probability density according to bins.
%   Returned array 'av' is the average observable array, at current time.
%   Return dimension is: [line index, 1, space index, (bin index),1].
%   The  bin index is used to index over probabilities, omitted otherwise.
%   Called by xpath,xpathfb,xpathw
%   xSPDE functions are licensed by Peter D. Drummond, (2020) - see License
%
%%%%%%  CHECK THE n-th OBSERVE FUNCTION AND MEASURE VARIABLES IF NOT EMPTY

if isempty (p.observe{n})                      %%check if observe defined
    av = {};                                   %%no observe defined
    return;                                    %%return if no function  
end                                            %%end observe function check
trans = p.transforms{n};                       %%store transform switch
if sum(trans(2:p.dimension)) > 0               %%if transform is needed
    a = xnft(a,trans,p);                       %%make a Fourier transform
end                                            %%end if transform needed
o = xshape(p.observe{n}(a,p),0,p);             %%Call the observe function
if p.thresholdw > 0 && p.d.obs{n}(end) > 1     %%If weighted averages
      Om = exp(a(p.fields,:));                 %%calculate weighting
      o = Om.*o;                               %%weight the averages
end                                            %%End if weighted averages

%%%%%%  GET SCATTER PLOT IF NEEDED

if p.scatters{n} > 0                           %%Check for scatter plot
    o = reshape(o,[p.d.obs{n}]);               %%Reshape for averaging
    av = zeros(p.scatters{n},p.nspace,1);      %%Initialize scatter plot
    for lines = 1:p.scatters{n}                %%Index over the lines
      for d1 = 1:p.d.obs{n}(1)
        av(lines+(d1-1)*p.scatters{n} ,:,1) = o(d1,:,lines);%%Store data
      end
    end                                        %%End index over the lines
    return;                                    %%End the xdata function
end                                            %%End scatter plot

%%%%%%   GET PROBABILITIES IF BINS AVAILABLE
                                             
if ~isempty(p.binranges{n})                    %%If probabilities needed
    binranges = p.binranges{n};                %%Store the n-th bins
    or  = real(o);                             %%Get real part of data
    sz  = size(or);                            %%Get size of data  
    k   = ones(1,sz(2));                       %%Initialize total indices
    in  = ones(1,sz(2));                       %%Initialize inrange switch
    tot = 1;                                   %%Initialize cumulative size
    for m = 1:length(binranges)                %%Loop over binned variables
        pl = p.d.bins{n}(m);                   %%Get m-th dimension
        b1 = binranges{m}(1);                  %%Get start of m-th bin
        db = binranges{m}(2) - binranges{m}(1);%%Get delta of m-th bin
        om = floor((or(m,:)- b1)/db);          %%Calculate m-th index  
        inrange = (om>=0).*(om<pl);            %%Is m-th index in range?
        om = om.*inrange;                      %%Make m-th index in range
        k  = k + om*tot;                       %%Compute total indices
        tot = tot*pl;                          %%Get cumulative bin size
        in = inrange.*in/db;                   %%All indices in range? 
    end                                        %%End loop over variables
    in=in/p.ensembles(1);                      %%Normalise probabilities
    indsize(1) = sz(2)/p.ensembles(1);         %%Index dimension 1
    av  = zeros(indsize(1),tot);               %%Reshape av for binning
    indsize(2) = p.ensembles(1);               %%Index dimension 2
    k = reshape(k,indsize);                    %%Reshape index array
    in = reshape(in,indsize);                  %%Reshape range switch array
    for s = 1:indsize(1)                       %%Loop over first index
      for e = 1:indsize(2)                     %%Loop over ensemble index
        av(s,k(s,e)) = av(s,k(s,e)) + in(s,e); %%Increment the probability
      end                                      %%End loop over ensemble
    end                                        %%End loop over first index
    return;                                    %%End the xdata function
end                                            %%End if probabilities
 
%%%%%% GET AVERAGE IF NEEDED

o = reshape(o,[p.d.obs{n}]);                   %%Reshape for averaging
if p.ensembles(1) > 1                          %%Check if ensembles
    av = mean(o,3);                            %%Take average if possible
    if p.thresholdw > 0 && p.d.obs{n}(3) > 1
      Om = reshape(Om,[1,1,p.d.obs{n}(3)]);    %%Reshape for averaging     
      av = av/mean(Om,3);                      %%Normalise
    end
else                                           %%No average is possible
    av = o;                                    %%Store current observable
end                                            %%End if ensembles
end                                            %%End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XDATA FUNCTION