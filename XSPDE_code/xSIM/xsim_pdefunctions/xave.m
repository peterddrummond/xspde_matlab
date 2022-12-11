function o  =  xave(o,varargin)               
%   a  =  XAVE(o[,av],r) averages a field over a spatial lattice. 
%   Input:   xave(o,r): variable 'o', parameters 'r'. 
%   Optional inputs: xave(o,av,r): averaging switch `av',
%   Output: lattice average of 'o'.
%   If 'r' only is present, averages only over the space-time lattice.
%   If 'av', 'r'  present, averages in directions where 'av(i)' > 0.
%   The av switch has a first dimension of time following xspde standard.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT ARGUMENTS

switch nargin                                    %%Check for arguments
    case 3                                       %%Case of three arguments
        av =varargin{1};                         %%av is input
        r = varargin{2};                         %%r is the last argument
    case 2                                       %%Case of two arguments
        r = varargin{1};                         %%r is the second argument
        av = ones(1,r.dimension);                %%dx is the volume element
    otherwise                                    %%Throw an error message
        error ('xSPDE error: xave requires 2 or 3 input arguments.')
end                                              %%End checks for arguments
o = xshape(o,r.indext,r);                        %%Reshape to standard
maxind =  numel(size(o))-r.indext;               %%Adjust input dimension
maxind = min([maxind,length(av),r.dimension]);   %%Get maximum dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  AVERAGE OVER DIMENSIONS 

for i = 2-r.indext:maxind                        %%Loop on good dimensions 
    if av(i)>0                                   %%If i-th average needed 
        o = mean(o,i+r.indext);                  %%Get average of o 
    end                                          %%End if i-th integration
end                                              %%End loop on dimension
end                                              %%End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XAVE FUNCTION