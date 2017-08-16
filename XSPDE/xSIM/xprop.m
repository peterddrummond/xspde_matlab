function a  =  xprop(a,r) 
%   a = XPROP(a,r) propagates a step in time for linear propagarion. 
%   Input: field a, lattice r.
%   Uses r.propagator to propagate in k-space
%   Imposes Dirichlet boundaries in x-space if needed.
%   Output: new field a. 
%   Unpacked first three dimensions of a are components, ensembles, time.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 
                                                            
if ~isequal(r.propagator,1)                   %%change required
  dmax = r.dimension+2;                       %%maximum dimension
  if ~isequal(r.propagator,0)                 %%k-space
    a =reshape(a,r.d.fields);                 %%reshape to array
    for nd = 4:dmax                           %%loop over dimension
        a = fft(a,[],nd);                     %%take Fourier transform
    end                                       %%end loop over dimension
    a = r.propagator.*a;                      %%propagate in Fourier space
    for nd = 4:dmax                           %%loop over dimension
        a = ifft(a,[],nd);                    %%inverse Fourier transform
    end                                       %%end loop over dimension
  else                                        %%x-space
    shape = [1,1,1];
    for nd = 4:dmax
        dir = nd-2;
        shape(1) =  prod(r.d.fields(1:nd-1)); 
        shape(2) =  r.d.fields(nd); 
        shape(3) =  prod(r.d.fields(1+nd:end)); 
        a = reshape(a,shape);                 %%Unflatten lattice
        if r.boundaries(dir,1) == 1           %%test for Dirichlet case
            a(:,1,:) = 0;                     %%vanishing boundaries
        end
        if r.boundaries(dir,2) == 1           %%test for Dirichlet case
            a(:,end,:) = 0;                   %%vanishing boundaries
        end
    end
  end                                         %%end if kspace
  a = reshape(a, r.d.a);                      %reshape to matrix
end                                           %%end if change
end                                           %%end xprop function