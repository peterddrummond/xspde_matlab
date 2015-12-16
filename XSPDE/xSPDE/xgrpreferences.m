function input = xgrpreferences (input)  
%   input  =  XGRPREFERENCES(input) sets default values for  graphics inputs.
%   Input:  'input' cell array.
%   Output: 'input' cell array with default graphics values. 
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

input = xinpreferences(input);                   %%check input
sequence = length(input);
for s = 1:sequence                               %%loop over sequence length
  in = input{s};                                 %%get input structure
  in.gversion =   xprefer(in,'gversion','1.04');
  in.minbar =     xcprefer(in,'minbar',in.graphs,{0.01});
  in.font =       xcprefer(in,'font',in.graphs,{18});
  in.headers =    xcprefer(in,'headers',in.graphs,{1});
  in.images =     xcprefer(in,'images',in.graphs,{0});  
  in.imagetype =  xcprefer(in,'imagetype',in.graphs,{1});
  in.transverse = xcprefer(in,'transverse',in.graphs,{0}); 
  in.pdimension = xcprefer(in,'pdimension',in.graphs,{4});
  in.xlabels =    xcprefer(in,'xlabels',4,{'t','x','y','z'});
  in.klabels =    xcprefer(in,'klabels',4,{'\omega','k_x','k_y','k_z'}); 
  in.compare =    xcprefer(in,'compare',in.graphs,{''});
  if in.print > 2
      in                                         %% show graph parameters
  end
  input{s} = in;                                 %%return cells
end
end                                              %% end function

%%Revision history:
%%0.73: 
%% added a function to test defaults, changed default labels, 
%% displays current graphics parameters
%% added a minimum size for error bars, minbar
%%1.0: 
%%Added hdf5file input

