function o  =  xint(o,varargin)                 
%   a  =  XINT(o,[dx,] r) integrates an input observable over a lattice. 
%   Input:  lattice variable 'o', averaging switch `dx', lattice 'r'.
%   Output: space integral of 'o' returned in all lattice points.
%   Integrates only in directions where `dx(i)' > 0.
%   First observable dimension is for ensembles, and is not integrated.
%   If no dx is present, then in.dx is assumed to apply.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

if nargin == 2
    r = varargin{1};                         %% r is the second argument
    dx = r.dx;                               %% dx is the volume element
else
    dx =varargin{1};                         %% dx is input
    r = varargin{2};                         %% r is the last argument
end
o = reshape(o,r.d.int);                      %%Unflatten lattice 
for i = 2:length(dx)
    if dx(i)>0
        index = ones(1,r.dimension);
        index(i) = r.d.int(i);
        o =repmat (sum(o,i)*dx(i),index);
    end
end
o  = reshape(o,r.d.r);                      %Flatten lattice 
end                                         %%end function
%  v1.05
%  Changed ones(1,4) to ones(1,r.dimension) for generality