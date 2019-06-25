function [a,boundvalue]  =  xsetbound(a,r) 
%   a = XSETBOUND(a,r) sets field boundary conditions for propagation. 
%   Input: field a, lattice r.
%   Imposes Dirichlet and/or Neuman boundaries in space if needed.
%   All boundary types and values can be set individually per field index.
%   They apply to any transverse dimension, and can change dynamically.
%   Boundary types are defined through r.boundaries{d}(i,b) for
%   space dimension (d=2,3..), field index (i=1,2..) and boundary b=(1,2) 
%   Here d is the transverse space dimension index, ie, not counting time
%   (a) r.boundaries{d}(i,j)  = 0 gives the default, periodic 
%   (b) r.boundaries{d}(i,j)  = -1 gives Neumann, prescribed derivative
%   (c) r.boundaries{d}(i,j)  = 1 gives Dirichlet, prescribed field
%   Boundary values are set on field i through boundvalue(i,..,j,..).
%   This array is dynamically returned by r.boundfun{dir}(a,r)
%   The index list of the boundary values matches the field index list.
%   Values are set sequentially by the dimension, and can be overwritten
%   Input: field a, lattice r.
%   Output: new field a. 
%   xSPDE functions are licensed by Peter D. Drummond, (2017) - see License
  
sz = size(a);                                 %%get the input size
boundvalue = r.boundfun(a,r);                 %%get the boundary values
shape = [sz(1),1,1,1];                        %%initialize the shape
    for dir = 2:r.dimension                   %%loop over space dimension
      bound=r.boundaries{dir};                %%store boundary switch
      if ~isequal(bound, zeros(size(bound)))  %%check if boundaries exist
        shape(2) =  prod(r.d.int(1:dir-1));   %%dimensions<dir 
        shape(3) =  r.d.int(dir);             %%dimension=dir 
        shape(4) =  prod(r.d.int(dir+1:end)); %%dimensions>dir 
        a = reshape(a,shape);                 %%Unflatten field to shape
        shape(3) =  2;                        %%set shape of boundary value
        boundval = reshape(boundvalue{dir},shape);
        for i=1:sz(1)                         %%loop over field index
          if bound(i,1)  == 1                 %%If lower Dirichlet
            a(i,:,1,:) = boundval(i,:,1,:);
          end
          if bound  == 1                      %%If upper Dirichlet
            a(i,:,end,:) = boundval(i,:,end,:);
          end
        end                                   %%end loop over field index
      end                                     %%end check if boundaries
    end                                       %%end loop over dimension
    a = reshape(a,sz);                        %%reshape to input size
end                                           %%end function