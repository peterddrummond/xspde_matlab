function o  =  xave(o,varargin)               
%   a  =  XAVE(o,varargin) averages an input observable over a lattice. 
%   Input: variable 'o' on a flattened lattice, original dimension r.d.int. 
%   Optional inputs: xave(o,av,r): averaging switch `av', parameters 'r'.
%   Output: lattice average of 'o' returned in all lattice points.
%   If neither 'av', 'r' are present all lattice indices are averaged.
%   If 'r' only is present, averages over the space-time lattice
%   If 'av', 'r'  present, averages in all directions where 'av(i)' > 0.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

sz = size(o);                              %%get original size to return
if nargin  == 1
    ave = mean (o,2);                      %%Take mean over second index
    o=repmat(ave,1,sz(2));                 %%Repeat over original lattice
else
    if nargin  == 2
        r = varargin{1};                   %%Parameters are last argument
        av = ones(1,r.dimension);
    else
        r = varargin{2};                   %%Parameters are last argument
        av = varargin{1};                  %%Switch is the second argument
    end
    o = reshape(o,r.d.int);                %%Unflatten lattice
    maxind = min(r.dimension,length(av));  %% get maximum dimension 
    for i = 2:maxind                       %%Loop over average switch
        if av(i)>0                         %%If average switch is true
            index = ones(1,r.dimension+1); %%Initialize index vector
            index(i+1) = r.d.int(i+1);     %%Set index to dimension size
            o =repmat (mean(o,i+1),index); %Repeat mean over lattice 
        end                                %%End if average switch is true
    end                                    %%End loop over average switch
end
o  = reshape(o,sz);                        %%Flatten lattice to input size
end                                        %%end function