function into  =  xinc(o,varargin)                 
%   into  =  XINC(o,varargin) integrates and compresses input observables. 
%   Input: matrix 'o' on a flattened lattice, original dimension r.d.int. 
%   Forms: xinc(o,r),xinc(o,dx,r): switch `dx', parameters 'r'.
%   Output: lattice full or partial average of 'o', compressed.
%   If neither 'dx', 'r' are present, second index summed, returns matrix.
%   If 'r' present,  space index is integrated, returns reduced array.
%   If 'dx', 'r'  present, integrates in space directions where 'dx(i)' > 0.
%   Using 'dx(1)>0' is ignored.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

switch nargin                                %% checks for arguments
    case 2                                   %% case of two arguments
        r = varargin{1};                     %% r is the second argument
        dx = r.dx;                           %% dx is the volume element
    case 3                                   %% case of three arguments
        dx =varargin{1};                     %% dx is input
        r = varargin{2};                     %% r is the last argument
    otherwise                                %% throw an error message
        error ('xSIM error: xinc requires 2 or 3 input arguments.')
end                                          %% end checks for arguments
maxind = min(r.dimension,length(dx));        %% get maximum dimension
s = size(o);                                 %% input size 
o = reshape(o,[s(1),r.d.int]);               %% unflatten lattice 
for i = 2:maxind                             %% loop to maximum dimension 
    if dx(i)>0                               %% is i-th integration needed? 
        into = sum(o,i+1)*dx(i);             %% integrate 
    end                                      %% end i-th integration
end                                          %% end loop to max dimension 

%  v1.2
%  Returns reduced array size