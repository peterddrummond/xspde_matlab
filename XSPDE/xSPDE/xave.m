function o  =  xave(o,varargin)               
%   a  =  XAVE(o,varargin) averages an input observable over a lattice. 
%   Input: variable 'o' on a flattened lattice, original dimension r.d.int. 
%   Optional inputs: xave(o,dx,r) with averaging switch `dx',lattice 'r'.
%   Output: lattice average of 'o' returned in all lattice points if no 'dx'.
%   If 'dx', 'r' not present even the ensemble direction is averaged.
%   If 'dx', 'r'  present, averages in all directions where 'dx(i)' > 0.
%   In the switched case, using 'dx(1)>0' averages over samples.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

len = length(o);
switch nargin
    case 1
        av = mean (o,2);
        o=repmat(av,len,1);
    case {2,3}
    r = varargin{2};
    o = reshape(o,r.d.int);                  %%Reshape to 4D sample+lattice
    dx = varargin{1};       
    for i = 1:length(dx)
        if dx(i)>0
            index = ones(1,4);
            index(i) = r.d.int(i);
            o =repmat (mean(o,i),index);
        end
    end
    o  = reshape(o,r.d.r);
end
end                                         %%end function