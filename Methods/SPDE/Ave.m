function o  =  Ave(o,varargin)               
%   a  =  AVE(o[,av],p) averages a field over a spatial lattice. 
%   Input:   Ave(o,p): variable 'o', parameters 'p'. 
%   Optional inputs: Ave(o,av,p): averaging switch `av',
%   Output: lattice average of 'o'.
%   If 'p' only is present, averages only over the space-time lattice.
%   If 'av', 'p'  present, averages in directions where 'av(i)' > 0.
%   The av switch has a first dimension of time following xspde standard.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT ARGUMENTS

switch nargin                                    %%Check for arguments
    case 3                                       %%Case of three arguments
        av =varargin{1};                         %%av is input
        p = varargin{2};                         %%r is the last argument
    case 2                                       %%Case of two arguments
        p = varargin{1};                         %%r is the second argument
        av = ones(1,p.dimensions);               %%dx is the volume element
    otherwise                                    %%Throw an error message
        error ('xSPDE error: xave requires 2 or 3 input arguments.')
end                                              %%End checks for arguments
maxind =  numel(size(o))-p.indext;               %%Adjust input dimension
maxind = min([maxind,length(av),p.dimensions]);  %%Get maximum dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  AVERAGE OVER DIMENSIONS

for i = 2-p.indext:maxind                        %%Loop on good dimensions 
    if av(i)>0                                   %%If i-th average needed 
        o = mean(o,i+p.indext);                  %%Get average of o 
    end                                          %%End if i-th integration
end                                              %%End loop on dimension
end                                              %%End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XAVE FUNCTION