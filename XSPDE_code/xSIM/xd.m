function d  =  xd(o,varargin) 
%   a = XD(o,D,r) calculates a spatial derivative of a field 
%   Uses periodic boundary conditions and Fourier transforms 
%   Input: field a, [derivative grid D,] lattice r.
%   Output: derivative  d. 
%   If no input derivative specified, takes derivative in x direction
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
                                                            
sz = size(o);                            %%get original size to return
if nargin == 2
    r = varargin{1};                     %% r is the second argument
    D = r.Dx;                            %% r.Dx is the derivative
else
    D =varargin{1};                      %% D is input
    r = varargin{2};                     %% r is the last argument
end
o = reshape(o,[sz(1),r.d.int]);          %%Unflatten lattice
for nd = 2:r.dimension                   %%loop over dimension
        o = fft(o,[],nd+1);              %%take Fourier transform
end                                      %%end loop over dimension
o = reshape(o, [sz(1),r.nlattice]);
for n  =  1: sz(1)
    o(n,:) = D.*o(n,:);                  %%derivative in Fourier space
end
o =reshape(o,[sz(1),r.d.int]);           %%reshape to array 
for nd =  2:r.dimension                  %%loop over dimension
    o = ifft(o,[],nd+1);                 %%inverse Fourier transform
end                                      %%end loop over dimension
d = reshape(o, sz );                     %%reshape to matrix
end                                      %%end xd function