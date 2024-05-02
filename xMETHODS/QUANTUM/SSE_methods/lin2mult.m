function [mult] = lin2mult(sz,lin)
%LIN2MULT Multiple subscript vector from a linear index.
len =    length(sz);
mult   = ones(1,len);
k      = cumprod(sz);
for i = len:-1:2
        vi = rem(lin-1, k(i-1)) + 1;
        mult(i) = (lin - vi)/k(i-1) + 1;
        lin = vi;
end
mult(1)=lin;
end