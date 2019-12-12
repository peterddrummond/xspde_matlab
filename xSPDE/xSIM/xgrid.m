function r = xgrid(n,r)                       
%   r = XGRID (r)  sets grid points in lattice from coordinate vectors.
%   Input:  struct r with coordinates.
%   n=1 for k-space in propagate order, n=2 for k-space in graph order
%   Output: struct r including all grid points in space and momentum. 
%   Note, for derivatives,  r.D{d} = 1i*k{d}.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

e = 1:r.ensembles(1);                    
d = r.dimension;
if n == 1
    r.r = cell(1,max(4,d));                         %%cell of position-grids
    r.k = cell(1,max(4,d));                         %%cell of k-grids
    r.D = cell(1,max(4,d));                         %%cell of derivatives
    [~,r.r{2:d}] = ndgrid(1,r.xc{2:d},e);           %%make x-grids
    [~,r.k{2:d}] = ndgrid(1,r.kcp{2:d},e);          %%make k-propagate grids
    for id = 2:d                                    %%loop over dimension
        r.r{id} = reshape (r.r{id},r.d.r);
        r.k{id} = reshape (r.k{id},r.d.r);
        r.D{id} = 1i*r.k{id};
    end 
    [r.x,r.y,r.z] =    deal(r.r{2:4});              %%cell of x-grids
    [r.kx,r.ky,r.kz] = deal(r.k{2:4});              %%cell of k-grids
    [r.Dx,r.Dy,r.Dz] = deal(r.D{2:4});              %%cell of k-grids  
    if r.print == 2
      display(r.x)
      display(r.kx);
    end
    %if (r.noises(2)+r.randoms(2)) > 0               %%if k-noise required
     % Kr =  r.rfilter(r);                       %%get k random filters
      %Kn =  r.nfilter(r);                       %%get k noise filters
      %r.infilt = reshape(Kr,[r.randoms(2),r.d.int]); %%reshape input filter
      %r.noisefilt = reshape(Kn,[r.noises(2),r.d.int]);%%reshape filter
    %end                                             %%end if k-noise 
else 
    for id = 2:d                            %%loop over dimension           
        r.k{id} = fftshift(r.k{id},1+d);
    end
    r.kx  =  fftshift(r.kx,3);
    if d>2
         r.ky  =  fftshift(r.ky,4);
         if d>3
              r.kz  =  fftshift(r.kz,5);
         end
   end
end
end