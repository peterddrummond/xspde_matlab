function a  =  xEuler(a,xi,r)       
%   a = XEULER(a,xi,dt,r)  propagates a step using the Euler method.   
%   Input:  field 'a',  noise 'xi', lattice 'r'.
%   Output: new field 'a'.  Output is reshaped by r.prop to flat matrix.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

     a = a+reshape(r.da(a,xi,r)*r.dtr,r.d.a);
     a = r.prop(a,r);                            %%Propagate linear term
end                                                      