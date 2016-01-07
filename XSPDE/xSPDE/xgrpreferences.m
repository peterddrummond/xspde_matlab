function input = xgrpreferences (input)  
%   input  =  XGRPREFERENCES(input) sets default values for  graphics inputs.
%   Input:  'input' cell array.
%   Output: 'input' cell array with default graphics values. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

input = xinpreferences(input);                   %%check input
sequence = length(input);
for s = 1:sequence                               %%loop over sequence length
  in = input{s};                                 %%get input structure
  in.gversion =   xprefer(in,'gversion','xGRAPH1.05');
  in.minbar =     xcprefer(in,'minbar',in.graphs,{0.01});
  in.font =       xcprefer(in,'font',in.graphs,{18});
  in.headers =    xcprefer(in,'headers',in.graphs,{1});
  in.images =     xcprefer(in,'images',in.graphs,{0});  
  in.imagetype =  xcprefer(in,'imagetype',in.graphs,{1});
  in.transverse = xcprefer(in,'transverse',in.graphs,{0}); 
  in.pdimension = xcprefer(in,'pdimension',in.graphs,{4});
  in.compare =    xcprefer(in,'compare',in.graphs,{''});
  if in.numberaxis
    in.xlabels =    xcprefer(in,'xlabels',in.dimension,{'t','x_1','x_2','x_3'});
    in.klabels =    xcprefer(in,'klabels',in.dimension,{'\omega','k_1','k_2','k_3'}); 
  else
    in.xlabels =    xcprefer(in,'xlabels',in.dimension,{'t','x','y','z'});
    in.klabels =    xcprefer(in,'klabels',in.dimension,{'\omega','k_x','k_y','k_z'}); 
  end
  if in.print > 2
      display(in,'graph parameters');            %% show graph parameters
  end
  input{s} = in;                                 %%return cells
end
end                                              %% end function


