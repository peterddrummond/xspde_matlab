function a  =  xEuler(a,xi,r)       
%   a = XEULER(a,xi,dt,r)  propagates a step using the Euler method.   
%   Input:  field 'a',  noise 'xi', lattice 'r'.
%   Output: new field 'a'.  Output is reshaped by r.prop to flat matrix.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

     a(1:r.fields,:) = r.prop(a(1:r.fields,:),r);  %%Propagate linear term
     a(1:r.fields,:) = a(1:r.fields,:)+reshape(r.da(a,xi,r)*r.dtr,r.d.a);
     if r.defines
        a = [a(1:r.fields,:);r.define(a,xi,r)];     %%Get defined fields
     end
end                                                      