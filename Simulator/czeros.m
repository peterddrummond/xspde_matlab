function o = czeros(a,varargin)
%   o = CZEROS (a(,b)) makes cells of zeros, size [a{c},b],
%   where a is a cell array giving sizes depending on c, 
%   b is a fixed size which is ignored if void or missing
%   xSPDE functions are licensed by Peter D. Drummond, (2023) - see License

len = length(a);                                 %% Initialize length
o   = cell(1,len);
b = [];
if nargin > 1
    b = varargin{1};
end
for c = 1:len
      o{c} = zeros([a{c},b]);                    %%Set fields as zeros
end
end                                              %% End function