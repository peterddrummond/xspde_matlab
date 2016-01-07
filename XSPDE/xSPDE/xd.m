function d  =  xd(ac,varargin) 
%   a = XD(ac,D,r) calculates a spatial derivative of a field component 
%   Input: field a, [derivative grid D,] lattice r.
%   Output: derivative  d. 
%   If no input derivative specified, takes derivative in x direction
%   Note, first dimension of ac is r.ensembles(1).
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 
                                                            
if nargin == 2
    r = varargin{1};                         %% r is the second argument
    D = r.D.x;                               %% r.D.x is the derivative
else
    D =varargin{1};                          %% D is input
    r = varargin{2};                         %% r is the last argument
end
ac =reshape(ac,r.d.int);                     %%reshape to array 
    for nd = 2:r.dimension                   %%loop over dimension
        ac = fft(ac,[],nd);                  %%take Fourier transform
    end                                      %%end loop over dimension
    ac = D.*ac;                              %%derivative in Fourier space
    for nd =  2:r.dimension                  %%loop over dimension
        ac = ifft(ac,[],nd);                 %%inverse Fourier transform
    end                                      %%end loop over dimension
    d = reshape(ac, [1,r.nlattice]);        %%reshape to matrix
end                                          %%end xd function