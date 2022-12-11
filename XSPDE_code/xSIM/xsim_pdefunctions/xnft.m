function a  =  xnft(a,varargin)                 
%   a  =  XNFT(o,[dx,] r) Fourier transforms a field with normalization. 
%   Input:  field  variable 'a', switch `trans', lattice 'r'.
%   First field index is not transformed.
%   In d space dimensions, 'trans' has length equal to r.dimension = d+1.
%   Output: space Fourier transform of 'a' returned in all lattice points.
%   Transforms in directions where `trans(i)' > 0, ignoring 'trans(1)'.
%   If no trans is present, Fourier transforms on all space of dimension d,
%   out(k) = 1/(2*pi)^(d/2)*\int dx (a(x)exp -ikx)
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2017) - see License

switch nargin                                %% checks for arguments
    case 3                                   %% case of three arguments
        trans =varargin{1};                  %% dx is input
        r = varargin{2};                     %% r is the last argument
    case 2                                   %% case of two arguments
        r = varargin{1};                     %% r is the second argument
        trans = ones(1,r.dimension);         %% trans is transform switch
    otherwise                                %% throw an error message
        error ('xSPDE error: xnft requires 2 or 3 input arguments.')
end                                          %% end checks for arguments
maxind = min(r.dimension,length(trans));     %% get maximum dimension
for nd = 2:maxind                            %% loop over space dimension 
    if trans(nd) > 0                         %%if FFT required
         a = fft(a,[],nd)*r.kfact(nd);       %%Fourier transform
    end                                      %%end if FFT required
end                                          %% end loop to max dimension 
end                                          %%end function