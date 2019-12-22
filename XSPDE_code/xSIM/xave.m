function o  =  xave(o,varargin)               
%   a  =  XAVE(o,varargin) averages a field over a spatial lattice. 
%   Input: variable 'o'; if matrix, second dimension is averaged. 
%   Optional inputs: xave(o,av,r): averaging switch `av', parameters 'r'.
%   Output: lattice average of 'o' returned in all lattice points.
%   If neither 'av', 'r' are present, last lattice index is averaged also.
%   If 'r' only is present, averages over the space-time lattice.
%   If 'av', 'r'  present, averages in directions where 'av(i)' > 0.
%   If av(NDIMS+1) > 0, averages over the  ensemble(1) members.
%   The av switch has first dimension of time following xspde standard.
%   Note that in most xSPDE observables, there is only one time-point.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT ARGUMENTS

sz = size(o);                                %% Get original size to return
switch nargin                                %% Check for arguments
    case 3                                   %% Case of three arguments
        av =varargin{1};                     %% av is input
        r = varargin{2};                     %% r is the last argument
    case 2                                   %% Case of two arguments
        r = varargin{1};                     %% r is the second argument
        av = ones(1,r.dimension);            %% dx is the volume element
    case 1                                   %% Case of one argument
        ave = mean (o,2);                    %% Take mean over second index
        o=repmat(ave,1,sz(2));               %% Repeat over  lattice
        return;                              %% Return to calling function
    otherwise                                %% Throw an error message
        error ('xSPDE error: xave requires 1 2 or 3 input arguments.')
end                                          %% End checks for arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RESHAPE DATA FOR AVERAGES  

if (length(sz) == 2)                         %% If standard matrix input
        o = reshape(o,[sz(1),r.d.int]);      %% Unflatten lattice
end                                          %% End if standard matrix
maxind = min(ndims(o)-1,length(av));         %% Get maximum dimension
sizeo = size(o);                             %% Get reshaped size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  AVERAGE OVER DIMENSIONS 

for i = 1:maxind                             %% Loop over all dimensions 
    if av(i)>0                               %% If i-th average needed 
        index = ones(1,ndims(o));            %% Initialize index vector
        index(i+1) = sizeo(i+1);             %% Set index to dimension size
        o1 = mean(o,i+1);                    %% Get average of o 
        o =repmat (o1,index);                %% Repeat mean over lattice 
    end                                      %% End if i-th integration
end                                          %% End loop to max dimension 
o  = reshape(o,sz);                          %% Restore original dimension
end                                          %% End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XAVE FUNCTION