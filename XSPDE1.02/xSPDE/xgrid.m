function r = xgrid (r)                       
%   r = XGRID (r)  sets grid points in lattice from coordinate vectors.
%   Input:  struct r with coordinates.
%   Output: struct r including all grid points. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

e = ones(1,r.ensembles(1));                     %%vector length = sample
r.ky =0;
r.kz =0;
switch r.dimension
  case 2
    [~,~,r.x] = ndgrid(1,e,r.xc{2});                           %%X grid
    r.x = reshape (r.x,r.d.r);
    [~,~,r.kx] = ndgrid(1,e,r.kc{2});                          %%K grid
    r.kx = reshape (r.kx,r.d.r);
  case 3
    [~,~,r.x,r.y] = ndgrid(1,e,r.xc{2},r.xc{3});               %%X grid
    r.x = reshape (r.x,r.d.r);
    r.y = reshape (r.y,r.d.r);
    [~,~,r.kx,r.ky] = ndgrid(1,e,r.kc{2},r.kc{3});             %%K grid
    r.kx = reshape (r.kx,r.d.r);
    r.ky = reshape (r.ky,r.d.r);
  case 4
    [~,~,r.x,r.y,r.z] = ndgrid(1,e,r.xc{2},r.xc{3},r.xc{4});   %%X grid      
    r.x = reshape (r.x,r.d.r);
    r.y = reshape (r.y,r.d.r);
    r.z = reshape (r.z,r.d.r);
    [~,~,r.kx,r.ky,r.kz] = ndgrid(1,e,r.kc{2},r.kc{3},r.kc{4});%%K grid
    r.kx = reshape (r.kx,r.d.r);
    r.ky = reshape (r.ky,r.d.r);
    r.kz = reshape (r.kz,r.d.r);
end
end
