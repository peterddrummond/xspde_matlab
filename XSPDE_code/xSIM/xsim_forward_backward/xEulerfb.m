function a  =  xEulerfb(a,pr,xi,r)       
%   a = XEULERFB(a,p,xi,r)  propagates a step using the Euler method,
%   as part of a Picard iteration solution of a forward-backward SDE.   
%   Input:  field 'a', previous iteration 'p', noise 'xi', lattice 'r'.
%   Previous iteration is time-reversed, and evaluated at future point!
%   Output: new field 'a'.  Output is reshaped by r.prop to flat matrix.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

     a = r.prop(a,r);               %%Propagate linear term
     p = pr{2};
     a= a + r.da(a,p,xi,r)*r.dtr;
end                                                      