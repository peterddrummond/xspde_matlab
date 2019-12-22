function a  =  xpropa(a,r) 
%   a = XPROPA(a,r) propagates a step in time for linear propagarion.
%   includes anti-aliasing.   
%   Input: field a, lattice r.
%   Uses r.propagator to propagate in k-space or x-space as necessary
%   Output: new field a. 
%   If r.propagator = 0 no action is taken
%   Unpacked first three dimensions of a are components, ensembles, time.
%   xSPDE functions are licensed by Peter D. Drummond, (2017) - see License 
                                                            
if r.propagator ~= 0                          %%No xprop required
  if r.dimension > 1                          %%PSDE
    if ~isequal(r.propagator,1)               %%If IP is needed
      a =reshape(a,r.d.fields);               %%reshape to array
      dmax = r.dimension+2;                   %%maximum dimension
      for nd = 3:dmax                         %%loop over dimension
        a = fft(a,[],nd);                     %%take Fourier transform
      end                                     %%end loop over dimension
      a = r.propagatora.*a;                   %% anti-aliasing propagation
      for nd = 3:dmax                         %%loop over dimension
        a = ifft(a,[],nd);                    %%inverse Fourier transform
      end                                     %%end loop over dimension
      a = reshape(a, r.d.a);                  %reshape to matrix
    end                                       %%end if IP
  else                                        %%SDE
    if ~isequal(r.propagator,1)               %%If IP is needed
      a = r.propagator.*a;                    %%propagate in x-space
    end
  end
end                                           %%end if r.propagator
if r.setboundaries
      [a,r.boundvalue]  =  xsetbound(a,r);   %%set boundary values
end
end                                           %%end xprop function