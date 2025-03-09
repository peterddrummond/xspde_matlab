function a  =  xnft(a,varargin)                 
%   a  =  XNFT(o,[dx,] r) Fourier transforms field cells with normalization 
%   Input:  field  cells 'a', switch `trans', lattice 'p'.
%   First field index is not transformed.
%   In d space dimensions, 'trans' has length equal to p.dimensions = d+1.
%   Output: space Fourier transform of 'a' returned in all lattice points.
%   Transforms in directions where `trans(i)' > 0, ignoring 'trans(1)'.
%   If no trans is present, Fourier transforms on all space of dimension d,
%   out(k) = 1/(2*pi)^(d/2)*\int dx (a(x)exp -ikx)
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2017) - see License

switch nargin                                %% checks for arguments
    case 3                                   %% case of three arguments
        trans = varargin{1};                 %% dx is input
        p = varargin{2};                     %% r is the last argument
    case 2                                   %% case of two arguments
        p = varargin{1};                     %% r is the second argument
        trans = ones(1,p.dimensions);        %% trans is transform switch
    otherwise                                %% throw an error message
        error ('xSPDE error: xnft requires 2 or 3 input arguments.')
end                                          %% end checks for arguments
maxind = min(p.dimensions,length(trans));    %% get maximum dimension
for c = 1:length(a)
  sz =  size(a{c});
  for nd = 2:maxind                          %% loop over space dimension 
    if trans(nd) > 0  && sz(nd) > 1          %%if FFT required
      a{c} = fft(a{c},[],nd)*p.kfact(nd);    %%Fourier transform
    end                                      %%end if FFT required
  end                                        %% end loop to max dimension
end
end                                          %%end function