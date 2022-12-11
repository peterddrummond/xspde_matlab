function r = xgrid(n,r)                       
%   r = XGRID (r)  sets grid points in lattice from coordinate vectors.
%   Input:  struct r with coordinates.
%   n=0 for k-space in propagate order
%   n>0 for x,k,o-space for graph (n)
%   Grid redundancy allows user choice of x, k, r labels for comparisons
%   Output: struct r including all grid points in space and momentum.
%   Assumes Matlab index broadcasting is available
%   Note, for derivatives,  r.D{d} = 1i*k{d}.
%   First dimension is the field index
%   Called by: xsim
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License

d = r.dimension;                                    %%space-time dimension
if n == 0                                           %%If n=0, packed grids
    r.r = cell(1,max(4,d));                         %%cell of x-grids
    r.k = cell(1,max(4,d));                         %%cell of k-grids
    r.D = cell(1,max(4,d));                         %%cell of derivatives
    [~,r.r{2:d}] = ndgrid(1,r.xc{2:d});             %%make x-grids
    for j = 1:r.fields                              %%loop over fields
      [~,k1{2:d}] = ndgrid(1,r.kcp{j}{2:d});        %%make shifted k-grids
      for id = 2:d                                  %%loop over dimension
        r.k{id}(j,:) = reshape(k1{id},[1,r.nspace]);%%reshape kgrid
      end                                           %%end loop on dimension
    end                                             %%end loop over fields
    for id = 2:d                                    %%loop over dimension
      r.r{id} = reshape (r.r{id},r.d.r);            %%pack lattice arrays
      bound = r.boundaries{id};                     %%store boundary switch
      if isequal(bound,zeros(size(bound)))
        r.k{id} = reshape(r.k{id}(1,:),[1,r.d.space]);   %%pack 1 k-array
      else
        r.k{id} = reshape(r.k{id},[r.fields,r.d.space]); %%pack all arrays
      end
      r.D{id} = 1i*r.k{id};                         %%packed derivative
    end                                             %%end loop on dimension
    [r.x,r.y,r.z] =    deal(r.r{2:4});              %%cell of x-grids
    [r.kx,r.ky,r.kz] = deal(r.k{2:4});              %%cell of k-grids
    [r.Dx,r.Dy,r.Dz] = deal(r.D{2:4});              %%cell of Dx-grids
else                                                %%If n>0,unpacked grids                                           %%If n>0, graphics
    xk = r.xk{n};                                   %%Initialize n-th axes 
    d = length(xk);                                 %%space-time dimension
    r.r = cell(1,max(4,d));                         %%cell of x-grids
    [~,r.r{1:d}] = ndgrid(1,xk{1:d});               %%make x-grids
    r.k = r.r;                                      %%Fourier space labels
    [r.t,r.x,r.y,r.z] =    deal(r.r{1:4});          %%cell of x-grids
    [r.w,r.kx,r.ky,r.kz] = deal(r.k{1:4});          %%cell of x-grids
    if r.print == 3                                 %%test if print needed
            display(r.x)                            %%print r-arrays
    end                                             %%end test if print
    
end                                                 %%end test n
end                                                 %%end grid function