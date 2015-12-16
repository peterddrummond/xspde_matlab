function a  =  xEuler(a,xi,dt,r)       
%   a = XEULER(a,xi,dt,r)  propagates a step using the Euler method.   
%   Input:  field 'a', lattice 'r', noise 'xi', step 'dt'.
%   Output: new field 'a'.  Output is reshaped by r.prop to flat matrix.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

a = a+reshape(r.da(a,xi,r)*dt,r.d.a);          %%derivative
a = r.prop(a,r);                               %%Propagate initial field
end   
                                                            
