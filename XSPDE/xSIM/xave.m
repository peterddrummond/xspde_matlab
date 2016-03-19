function o  =  xave(o,varargin)               
%   a  =  XAVE(o,varargin) averages an input observable over a lattice. 
%   Input: variable 'o' on a flattened lattice, original dimension r.d.int. 
%   Optional inputs: xave(o,av,r): averaging switch `av', parameters 'r'.
%   Output: lattice average of 'o' returned in all lattice points.
%   If neither 'av', 'r' are present all lattice indices are averaged.
%   If 'r' only is present, averages only over the local ensemble
%   If 'av', 'r'  present, averages in all directions where 'av(i)' > 0.
%   Using 'av(1)>0' averages over samples in the local ensemble.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

switch nargin                              %%Switch  on argument number
    case 1                                 %%Average over entire lattice
        len = length(o);                   %%Store the input length
        ave = mean (o,2);                  %%Take mean over second index
        o=repmat(ave,len,1);               %%Repeat over original lattice
    case 2                                 %%Average over local ensemble
        r = varargin{1};                   %%Parameters r = second argument
        o = reshape(o,r.d.int);            %%Unflatten lattice
        index = ones(1,r.dimension);       %%Initialize index vector
        index(1) = r.d.int(1);             %%Set first index to ensemble 
        o =repmat (mean(o,1),index);       %%Repeat mean over lattice 
        o  = reshape(o,r.d.r);             %%Flatten lattice
    case 3                                 %%Average in switch directions
        r = varargin{2};                   %%Parameters are last argument
        o = reshape(o,r.d.int);            %%Unflatten lattice
        av = varargin{1};                  %%Switch is the second argument
        maxind = min(r.dimension,length(av));%% get maximum dimension 
        for i = 1:maxind                   %%Loop over average switch
          if av(i)>0                       %%If average switch is true
            index = ones(1,r.dimension);   %%Initialize index vector
            index(i) = r.d.int(i);         %%Set index to dimension size
            o =repmat (mean(o,i),index);   %Repeat mean over lattice 
          end                              %%End if average switch is true
        end                                %%End loop over average switch
        o  = reshape(o,r.d.r);             %%Flatten lattice
    otherwise                              %% throw an error message
        error ('xSPDE error: xave has more than three arguments.')
end                                        %%End switch depending on number
end                                        %%end function
%  v1.06
%  Adds error message