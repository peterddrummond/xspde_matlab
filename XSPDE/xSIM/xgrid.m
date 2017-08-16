function r = xgrid(r)                       
%   r = XGRID (r)  sets grid points in lattice from coordinate vectors.
%   Input:  struct r with coordinates.
%   Output: struct r including all grid points in space and momentum. 
%   Note, for derivatives,  r.D{d} = 1i*k{d}.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

e = ones(1,r.ensembles(1));                    
d = r.dimension;
x = cell(1,4);                                  %%cell of x-grids
kc = cell(1,4);                                 %%cell of k-grids
kp = cell(1,4);                                 %%cell of k-grids
[~,~,x{2:d}] = ndgrid(e,1,r.xc{2:d});           %%make x-grids
[~,~,kc{2:d}] = ndgrid(e,1,r.kc{2:d});          %%make k-graphics-grids
[~,~,kp{2:d}] = ndgrid(e,1,r.kp{2:d});          %%make k-propagate grids
for id = 2:d                                    %%loop over dimension
        x{id} = reshape (x{id},r.d.r);
        kc{id} = reshape (kc{id},r.d.r);
        kp{id} = reshape (kp{id},r.d.r);
end 
if r.numberaxis||d>4
    r.x = cell(1,d);                            %%cell of x-grids
    r.k = cell(1,d);                            %%cell of k-grids   
    for id = 2:d                                %%loop over dimension
        r.x{id} = x{id};
        r.k{id} = kp{id};
    end 
else
    [r.x,r.y,r.z] =    deal(x{2:4});            %%cell of x-grids
    [r.kx,r.ky,r.kz] = deal(kp{2:4});           %%cell of k-propagate-grids   

end

if (r.noises(2)+r.randoms(2)) > 0               %%if k-noise required
      Kr =  r.rfilter(r);                       %%get k random filters
      Kn =  r.nfilter(r);                       %%get k noise filters
      r.infilt = reshape(Kr,[r.randoms(2),r.d.int]); %%reshape input filter
      r.noisefilt = reshape(Kn,[r.noises(2),r.d.int]);%%reshape filter
end                                             %%end if k-noise 

if r.numberaxis||d>4
    r.D = cell(1,d);                            %%cell of k-grids
    for id = 2:d                                %%loop over dimension
        r.D{id} = 1i*kp{id};
        r.k{id} = kc{id};
    end 
else
    [r.kx,r.ky,r.kz] = deal(kc{2:4});
    [r.Dx,r.Dy,r.Dz] = deal(kp{2:4});
    r.Dx  =  1i* r.Dx;
    r.Dy  =  1i* r.Dy;
    r.Dz  =  1i* r.Dz;
end
end