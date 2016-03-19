function r = xgrid(r)                       
%   r = XGRID (r)  sets grid points in lattice from coordinate vectors.
%   Input:  struct r with coordinates.
%   Output: struct r including all grid points in space and momentum. 
%   Note, for derivatives,  D{d} = 1i*k{d}.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

e = ones(1,r.ensembles(1));                         %%length = samples
d = r.dimension;
if r.numberaxis||d>4
    r.x = cell(1,d-1);                                %%cell of x-grids
    r.k = cell(1,d-1);                                %%cell of k-grids
    [~,~,r.x{1:d-1}] = ndgrid(1,e,r.xc{2:d});         %%make x-grids
    [~,~,r.k{1:d-1}] = ndgrid(1,e,r.kc{2:d});         %%make k-grids
    for id = 1:d-1                                    %%loop over dimension
        r.x{id} = reshape (r.x{id},r.d.r);
        r.k{id} = reshape (r.k{id},r.d.r);
    end 
else
  if d > 1
    x = cell(2,4);                                  %%cell of x-grids
    k = cell(2,4);                                  %%cell of k-grids
    [~,~,x{2:d}] = ndgrid(1,e,r.xc{2:d});           %%make x-grids
    [~,~,k{2:d}] = ndgrid(1,e,r.kc{2:d});           %%make k-grids
    for id = 2:d                                    %%loop over dimension
        x{id} = reshape (x{id},r.d.r);
        k{id} = reshape (k{id},r.d.r);
    end
    [r.x,r.y,r.z] =    deal(x{2:4});
    [r.kx,r.ky,r.kz] = deal(k{2:4});
  end
end
end
