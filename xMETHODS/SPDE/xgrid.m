function p = xgrid(n,p)                       
%   p = XGRID (p)  sets grid points in lattice from coordinate vectors.
%   Input:  struct p with coordinates, output index n.
%   n = 0 for k-space in propagate order
%   n > 0 for x,k,o-space for graph (n) and comparisons
%   Grid redundancy allows user choice of x, k, rp labels for comparisons
%   Output: struct p including all grid points in space and momentum.
%   Assumes Matlab index broadcasting is available
%   Note, for derivatives,  p.D{d} = 1i*k{d}.
%   First dimension is the field index
%   Called by: xsim, xpreferences
%   Licensed by Peter D. Drummond, (2024) - see License

if p.dimensions == 1 && (n==0 || p.bins{n} == 0)
    return
end
d = p.dimensions;                                %%space-time dimension
if n == 0                                        %%If n = 0, initial grids
  p.r = cell(1,max(4,d));                        %%cell of x-grids
  p.k = cell(p.fieldcells,max(4,d));             %%cell of k-grids
  [~,p.r{2:d}] = ndgrid(1,p.xc{2:d});            %%make x-grids
  for c=1:1 %p.fieldcells                           %%loop over cells
     p.d.space{c} = p.points{c}(2:d);            %%space points
     nspace = prod(p.d.space{c});                %%number of space points
    for j = 1:p.fields{c}                        %%loop over fields
      [~,k1{2:d}] = ndgrid(1,p.kcp{c,j}{2:d});   %%make shifted k-grids
      for id = 2:d                               %%loop over dimension
        p.k{c,id}(j,:) = reshape(k1{id},[1,nspace]);%%reshape kgrid
      end                                        %%end loop on dimension
    end                                          %%end loop over fields
  end                                            %%end loop over cells
  for id = 2:d                                   %%loop over dimension
    p.r{id} = reshape (p.r{id},p.d.r);           %%pack lattice arrays
    for c=1:1 %p.fieldcells                         %%loop over cells
      bound = p.boundaries{c,id};                %%store boundary switch
      if isequal(bound,zeros(size(bound)))
        p.k{c,id} = reshape(p.k{c,id}(1,:),[1,p.d.space{c}]);   
      else
        p.k{c,id} = reshape(p.k{id},[p.fields{c},p.d.space{c}]); 
      end
    end                                          %%end loop over cells
  end                                            %%end loop on dimension
  [p.x,p.y,p.z] =    deal(p.r{2:4});             %%cell of x-grids
  [p.kx,p.ky,p.kz] = deal(p.k{2:4});             %%cell of k-grids
else                                             %%If n>0,output grids                                           %%If n>0, graphics
  xk = p.xk{n};                                  %%Initialize n-th axes 
  d = length(xk);                                %%space-time dimension
  p.r = cell(1,max(4,d));                        %%cell of x-grids
  [~,p.r{1:d}] = xdgrid(1,xk{1:d});              %%make x-grids
  p.k = p.r;                                     %%Fourier space labels
  [p.t,p.x,p.y,p.z]    = deal(p.r{1:4});         %%cell of x-grids
  [p.w,p.kx,p.ky,p.kz] = deal(p.k{1:4});         %%cell of x-grids
  if p.verbose == 3                              %%test if print needed
      display(p.x)                               %%print r-arrays
  end                                            %%end test if print  
end                                              %%end test n
end                                              %%end grid function

function [varargout] = xdgrid(varargin)
%   p = XDGRID (n,p)  makes a broadcast lattice from coordinate vectors.
%   This is functionally similar to ndgrid, but reduces memory use
%   Only the n-th lattice index is filled for the n-th grid output
%   All other lattice indices are set to one, allowing broadcasting
%   xSPDE functions are licensed by Peter D. Drummond, (2024) - see License

d = nargout;                                     %%store number of outputs
varargout = cell(1,d);                           %%make output cells
for  j = 1:d                                     %%loop over output cells
    m = min(j,nargin);                           %%largest available input
    size = [ones([1,j-1]),length(varargin{m}),ones([1,d-j])];
    varargout{j} = reshape(varargin{m},size);
end                                              %%end loop over outputs
end                                              %%end xdgrid function