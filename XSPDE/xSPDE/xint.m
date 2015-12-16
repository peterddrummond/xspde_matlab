function o  =  xint(o,dx,r)                 
%   a  =  XINT(o,r,dx) integrates an input observable over a lattice. 
%   Input:  lattice variable 'o', lattice 'r', averaging switch `dx'.
%   Output: space integral of 'o' returned in all lattice points.
%   Integrates only in directions where `dx(i)' > 0.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

o = reshape(o,r.d.int);                    %%Reshape to 4D sample+lattice
for i = 2:length(dx)
    if dx(i)>0
        index = ones(1,4);
        index(i) = r.d.int(i);
        o =repmat (sum(o,i)*dx(i),index);
    end
end
o  = reshape(o,r.d.r);
end                                         %%end function
