function av = xdata(a,n,p)  
%   av = XDATA(a,n,p) stores data averages, scatters or probabilities
%   Input is the 'a' cell array of stochastic variables, stored in parallel
%   with dimensions of: [internal indices, space indices, ensemble index].
%   Input index n gives the number of the observe function.
%   Input struct p gives parameters and observe function handles.
%   Observe dimension is [line index, space indices, ensemble index].
%   If s = p.scatters{n} > 1, makes  scatter-plots of first s trajectories.
%   If p.bin{n} is not empty, stores probability density in bins.
%   Returned array 'av' is the average observable array, at current time.
%   Returned dimension is: [line index, 1, (space index)*(bin index)].
%   The  bin index is used to index over probabilities, or else omitted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by path
%   Calls     observe{n}
%   xSPDE functions licensed by Peter D. Drummond, (2024) - see License
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CHECK THE n-th OBSERVE FUNCTION AND MEASURE VARIABLES IF NOT EMPTY

if isempty (p.observe{n})                        %%check if observe defined
    av = [];                                     %%no observe defined
    return;                                      %%return if no function  
end                                              %%end function check
trans = p.transforms{n} ;                        %%store transform switch
if sum(trans(2:p.dimensions)) > 0                %%if transform is needed
    a = xnft(a,trans,p);                         %%make a Fourier transform
end                                              %%end if transform needed
obs = p.observe{n}(a{:},p);                      %%Call observe function
%obs = xshape(p.observe{n}(a{:},p),0,p);         %%Call the observe function
sz = size(obs);
obs = reshape(obs,[p.d.obs{n}]);                 %%Reshape for averaging

%%%%%%  GET SCATTER PLOT IF NEEDED

if p.scatters{n} > 1                             %%Check for scatter plot
    av = zeros(p.scatters{n},p.nspace,1);        %%Initialize scatter plot
    for lines = 1:p.scatters{n}                  %%Index over the lines
      for d1 = 1:p.d.obs{n}(1)                   %%Loop on obs dimenension
        av(lines+(d1-1)*p.scatters{n} ,:,1) = obs(d1,:,lines);%%Store data
      end                                        %%End loop on dimenension
    end                                          %%End index over the lines
    return;                                      %%Return from function
end                                              %%End if scatter plot

%%%%%%   GET PROBABILITIES IF BINS AVAILABLE
                                             
if ~isempty(p.binranges{n})                      %%If probabilities needed
     obs = reshape(obs,sz);                      %%Reshape for binning
    binranges = p.binranges{n};                  %%Store the n-th bins
    or  = real(obs);                             %%Get real part of data
    sz  = size(or);                              %%Get size of data  
    k   = ones(1,sz(2));                         %%Initialize total indices
    in  = ones(1,sz(2));                         %%Initialize inrange 
    tot = 1;                                     %%Initialize total size
    for m = 1:length(binranges)                  %%Loop over bins
        pl = p.d.bins{n}(m);                     %%Get m-th dimension
        b1 = binranges{m}(1);                    %%Get start of m-th bin
        db = binranges{m}(2) - binranges{m}(1);  %%Get delta of m-th bin
        om = floor((or(m,:)- b1)/db);            %%Calculate m-th index  
        inrange = (om>=0).*(om<pl);              %%Is m-th index in range?
        om = om.*inrange;                        %%Make m-th index in range
        k  = k + om*tot;                         %%Compute total indices
        tot = tot*pl;                            %%Get cumulative bin size
        in = inrange.*in/db;                     %%All indices in range? 
    end                                          %%End loop over variables
    in=in/p.ensembles(1);                        %%Normalise probabilities
    indsize(1) = sz(2)/p.ensembles(1);           %%Index dimension 1
    av  = zeros(indsize(1),tot);                 %%Reshape av for binning
    indsize(2) = p.ensembles(1);                 %%Index dimension 2
    k = reshape(k,indsize);                      %%Reshape index array
    in = reshape(in,indsize);                    %%Reshape range switch 
    for s = 1:indsize(1)                         %%Loop over first index
      for e = 1:indsize(2)                       %%Loop over ensemble index
        av(s,k(s,e)) = av(s,k(s,e)) + in(s,e);   %%Increment probability
      end                                        %%End loop over ensemble
    end                                          %%End loop over 1st index
else                                             %%Else not a probability
    av = mean(obs,4);                            %%Take average of ensemble
end                                              %%End if probabilities
end                                              %%End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XDATA FUNCTION