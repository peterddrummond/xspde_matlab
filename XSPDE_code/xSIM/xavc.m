 function o  =  xavc(o,varargin)               
%   o  =  XAVC(o,varargin) averages and compresses input observables. 
%   Input: matrix 'o' on a flattened lattice, original dimension r.d.int. 
%   Forms: xavc(o),xavc(o,r),xavc(o,av,r): switch `av', parameters 'r'.
%   Output: lattice full or partial average of 'o', compressed.
%   If neither 'av', 'r' are present, second index is averaged.
%   If 'r' present, unflattens, space index is averaged.
%   If 'av', 'r'  present, averages in all directions where 'av(i)' > 0.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

s = size(o);                               %% input size

switch nargin                              %%Switch  on argument number
    case 1                                 %%Average over entire lattice
        o = mean (o,2);                    %%Take mean over second index
    case 2                                 %%Average over all space
        r = varargin{1};                   %%Parameters r = second argument
        o = reshape(o,[s(1),r.d.int]);     %% unflatten lattice
        for i = 2:r.dimension              %%Loop over dimension
            o = mean(o,1+i);               %%Repeat mean over lattice
        end                                %%End loop 
    case 3                                 %%Average in switch directions
        r = varargin{2};                   %%Parameters are last argument
        o = reshape(o,[s(1),r.d.int]);     %% unflatten lattice
        av = varargin{1};                  %%Switch is the second argument
        maxind = min(r.dimension,length(av));%% get maximum dimension 
        for i = 2:maxind                   %%Loop over average switch
          if av(i)>0                       %%If average switch is true
            o = mean(o,1+i);               %%Repeat mean over lattice
          end                              %%End if average switch is true
        end                                %%End loop over average switch
    otherwise                              %% throw an error message
        error ('xSPDE error: xavc has more than three arguments')
end                                        %%End switch depending on number
end                                        %%end function
%  v1.2
%  Returns reduced array size
