function o  =  xave(o,varargin)               
%   a  =  XAVE(o,varargin) averages an input observable over a lattice. 
%   Input: variable 'o' on a flattened lattice, original dimension r.d.int. 
%   Optional inputs: xave(o,av,r) with averaging switch `av', and parameters 'r'.
%   Output: lattice average of 'o' returned in all lattice points.
%   If neither 'av', 'r' are present all lattice indices are averaged.
%   If 'av', 'r'  present, averages in all directions where 'av(i)' > 0.
%   In this case, using 'av(1)>0' averages over samples in the local ensemble.
%   If 'r' only is present, averages only over the local ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

len = length(o);                               %%Store the input bevtor length
switch nargin                                  %%Switch depending on argument number
    case 1                                     %%Average entire local lattice
        ave = mean (o,2);                      %%Take mean over second index
        o=repmat(ave,len,1);                   %%Repeat over original lattice
    case 2                                     %%Average only over local ensemble
        r = varargin{1};                       %%Parameters r are the second argument
        o = reshape(o,r.d.int);                %%Unflatten lattice to full dimension
        index = ones(1,r.dimension);           %%Initialize index vector
        index(1) = r.d.int(1);                 %%Set first index to ensemble size
        o =repmat (mean(o,1),index);           %%Repeat mean over original lattice 
        o  = reshape(o,r.d.r);                 %%Flatten lattice to internal dimension
    case 3                                     %%Average in average switch directions
        r = varargin{2};                       %%Parameters r are the last argument
        o = reshape(o,r.d.int);                %%Unflatten lattice to full dimension
        av = varargin{1};                      %%Average switch is the second argument      
        for i = 1:length(av)                   %%Loop over average switch
          if av(i)>0                           %%If average switch is true
            index = ones(1,r.dimension);       %%Initialize index vector
            index(i) = r.d.int(i);             %%Set index to dimension size
            o =repmat (mean(o,i),index);       %Repeat mean over original lattice 
          end                                  %%End if average switch is true
        end                                    %%End loop over average switch
        o  = reshape(o,r.d.r);                 %%Flatten lattice to internal dimension
end                                            %%End switch depending on argument number
end                                            %%end function