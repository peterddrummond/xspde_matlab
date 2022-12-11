function a  =  xprop(a,r) 
%   a = XPROP(a,r) propagates a step in time for linear propagation. 
%   Input: field a, parameters r.
%   Uses r.propagator to propagate in k-space if necessary
%   Output: new field a. 
%   If r.propagator = 0 no action is taken, apart from setting boundaries
%   xSPDE licensed by P. D. Drummond & S. Kiesewetter, (2022) - see License 
                                                            

if r.propagator ~= 0                             %%If propagator needed
  if r.dimension > 1                             %%If SPDE case
      sz = size(a);                              %%field size
      sz2 = prod(sz(2:end));                     %%slice size
      a = reshape(a,[sz(1),sz2]);                %%matrix field shape
      for nd = 2:r.dimension                     %%loop over dimension
        b = r.boundaries{nd};                    %%store boundary switch
        type  = b(:,1)+2*b(:,2);                 %%get integer switch
        for i=1:sz(1)                            %%loop over field index
          a1 = reshape(a(i,:),[1,sz(2:end)]);    %%take one slice
          %ab{nd,i,1} = a1(:,1,:);               %%future use
          %ab{nd,i,2} = a1(:,end,:);             %%future use
          switch type(i)                         %%switch on boundary type
            case -3                              %%Robin boundaries
              a1 = xdct1(a1,nd);                 %%take cosine transform
            case -1                              %%Dirichlet-Robin boundaries
              a1 = xdst3(a1,nd);                 %%take cosine transform
            case 0                               %%periodic boundaries
             a1 = fft(a1,[],nd);                %%take Fourier transform
             %a1 = xfft(a1,nd);                %%take Fourier transform
            case 1                               %%Robin-Dirichlet-boundaries
              a1 = xdct3(a1,nd);                 %%take cosine transform
            case 3                               %%Dirichlet boundaries
              a1 = xdst1(a1,nd);                 %%take full sine transform
            otherwise                            %%invalid boundary type
              error('Invalid IP boundary: [%d,%d]',b(i,1),b(i,2));
          end                                    %%end switch
          a(i,:) = reshape(a1,[1,sz2]);          %%restore one slice
        end                                      %%end loop over field
      end                                        %%end loop over dimension
      a = r.propagator.*reshape(a,sz);           %%propagate in k-space
      a = reshape(a,[sz(1),sz2]);                %%matrix field shape
      for nd = 2:r.dimension                     %%loop over dimension
        b = r.boundaries{nd};                    %%store boundary switch
        type  = b(:,1)+2*b(:,2);                 %%get integer switch
        for i=1:sz(1)                            %%loop over field index
          a1 = reshape(a(i,:),[1,sz(2:end)]);    %%take one slice
          switch type(i)                         %%switch on boundary type
            case -3                              %%Robin boundaries
              a1 = xdct1(a1,nd);                 %%take inverse transform
            case -1                              %%Dirichlet-Robin boundaries
              a1 = xdst2(a1,nd);                 %%take inverse transform
            case 0                               %%periodic boundaries
              a1 = ifft(a1,[],nd);               %%take Fourier transform
              %a1 = xifft(a1,nd);               %%take Fourier transform
            case 1                               %%Robin-Dirichlet boundaries
              a1 = xdct2(a1,nd);                 %%take inverse transform
            case 3                               %%Dirichlet boundaries        
              a1 = xdst1(a1,nd);                 %%take inverse transform
          end                                    %%end switch
          a(i,:) = reshape(a1,[1,sz2]);          %%restore one slice
        end                                      %%end loop over components
      end                                        %%end loop over dimension
      a = reshape(a,sz);                         %%restore the field shape
  else                                           %%else it is an SDE case
      a = r.propagator*a;                        %%propagate in time
  end                                            %%End if PSDE/SDE case 
end                                              %%end if propagator needed
if r.setboundaries                               %%If boundaries needed
   a  =  xsetbound(a,r);                         %%set boundary values
end                                              %%End if boundaries needed
end                                              %%end xprop function