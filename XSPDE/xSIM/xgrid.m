function r = xgrid(r)                       
%   r = XGRID (r)  sets grid points in lattice from coordinate vectors.
%   Input:  struct r with coordinates.
%   Output: struct r including all grid points in space and momentum. 
%   Note, for derivatives,  r.D{d} = 1i*k{d}.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

e = ones(1,r.ensembles(1));                    
d = r.dimension;
x = cell(1,4);                                  %%cell of x-grids
k = cell(1,4);                                  %%cell of k-grids
[~,~,x{2:d}] = ndgrid(e,1,r.xc{2:d});           %%make x-grids
[~,~,k{2:d}] = ndgrid(e,1,r.kc{2:d});           %%make k-grids
for id = 2:d                                    %%loop over dimension
        x{id} = reshape (x{id},r.d.r);
        k{id} = reshape (k{id},r.d.r);
end 
if r.numberaxis||d>4
    r.x = cell(1,d);                            %%cell of x-grids
    r.k = cell(1,d);                            %%cell of k-grids
    r.D = cell(1,d);                            %%cell of D-grids   
    for id = 2:d                                %%loop over dimension
        r.x{id} = x{id};
        r.k{id} = k{id};
        r.D{id} = 1i*r.k{id};
    end 
else
    [r.x,r.y,r.z] =    deal(x{2:4});
    [r.kx,r.ky,r.kz] = deal(k{2:4});
    r.Dx  =  1i* r.kx;
    r.Dy  =  1i* r.ky;
    r.Dz  =  1i* r.kz;
end
end
