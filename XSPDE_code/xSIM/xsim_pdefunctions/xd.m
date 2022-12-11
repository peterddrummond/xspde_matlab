function d  =  xd(o,varargin) 
%   a = XD(o,[D,]r) calculates a spatial derivative of a field 
%   Uses periodic boundary conditions and Fourier transforms 
%   Input: field a, [derivative grid D,], lattice r.
%   Output: derivative  d. 
%   If no input derivative specified, takes derivative in x direction
%   Contextual: if r.indext = 1 can differentiate over time, otherwise not
%   First dimension is the field index or line index
%   xSPDE functions are licensed by Peter D. Drummond, (2022) - see License
                                                            
if nargin == 2
    r = varargin{1};                         %% r is the second argument
    D = r.Dx;                                %% r.Dx is the derivative
else
    D = varargin{1};                         %% D is input
    r = varargin{2};                         %% r is the last argument
end
o = xshape(o,r.indext,r);                    %%Unflatten lattice if needed
if r.indext                                  %%If time index present
    D = reshape(D,[1,size(D)]);              %%Reshape differential
end                                          %%End if time index present
for nd = (2-r.indext):r.dimension            %%loop over dimension
        o = fft(o,[],nd+r.indext);           %%take Fourier transform
end                                          %%end loop over dimension
d = D.*o;                                    %%derivative in Fourier space
for nd =  (2-r.indext):r.dimension           %%loop over dimension
    d = ifft(d,[],nd+r.indext);              %%inverse Fourier transform
end                                          %%end loop over dimension
end                                          %%end xd function