    function input = xgpreferences (input,oldinput)  
%   input  =  XGPREFERENCES(input) sets default values for  graphics inputs.
%   Input:  'input' cell array, and any previous input 'oldinput'.
%   Output: 'input' cell array with updated and default graphics values.
%   Called by: xgraph
%   Needs:     xprefer, xcprefer
%   xGRAPH functions are licensed by Peter D. Drummond, (2021) - see License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  UPDATE OLD SEQUENCES
%                          
if ~isempty(oldinput)                            %% Check  oldinput ~empty
    oldsequence = length(oldinput);              %%get old sequence length
    for s= 1:oldsequence                         %%loop over old sequence 
        in = oldinput{s};                        %%get old sequence input
        fname = fieldnames(in);                  %%get old sequence  labels
        if s<= length(input)                     %%check new input length
          for j = 1:length(fname)                %%loop over old labels
            label = fname{j} ;                   %%set new label = old 
            if ~isfield(input{s},label)          %% If no label data 
              input{s}.(label) = in.(label);     %% Set new input to old
            end                                  %% End if no label data
          end                                    %% End for loop
        else                                     %%check no new input
          input{s} = oldinput{s};                %%set new input to old
        end                                      %% End s < = length
    end                                          %% End sequence loop
end                                              %% End if oldinput ~empty
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET THE STANDARD DEFAULT VALUE
%                          
sequence = length(input);                        %%get sequence length
linestyle = {'-k','--k','-.k',':k','-ok','--ok','-.ok',':ok','-+k','--+k'};
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE WHOLE SEQUENCE
%
%
for s = 1:sequence                               %%loop over sequence 
  in = input{s};                                 %%current inputs
  dgraphs = in.dgraphs; 
  octave = exist('OCTAVE_VERSION','builtin');    %%finds if octave version
  if ~isfield(in,'dimension')                    %% If no label data
      in.dimension = 0;
  end
  if ~isfield(in,'numberaxis')                   %% If no label data
      in.numberaxis = 0;
  end
  if ~isfield(in,'functions')                    %% If no functions data
      in.functions = 0;
  end
  if ~isfield(in,'graphs')                       %% If no graphs specified
      in.graphs = 1:dgraphs;                     %% Set to default value
  end
  maxd = max(in.dimension+1,in.maxd);            %%max number of dimensions
  dgraphs = max(in.graphs);                      %%max graph number
  xlabels=cell(1,in.dimension);                  %%cell array for xlabels 
  klabels=cell(1,in.dimension);                  %%cell array for klabels
  bounds{1} = num2cell(zeros(1,in.dimension));   %%initialize bounds
  if in.numberaxis || in.dimension > 4           %%set the default labels
      for i=1:in.dimension
          xlabels{i} = sprintf('r_%d',i);
          klabels{i} = sprintf('k_%d',i);
       end 
  else
      if in.dimension
        xlabels={'t','x','y','z'};
        klabels={'\omega','k_x','k_y','k_z'};
      else
          for nd = 1:maxd                   %% Loop over dimension
              label = sprintf('m_%d',nd);
              xlabels{nd}=label;
          end
      end
  end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET COMMON GRAPHICS DEFAULTS
%                          
  in.octave =     xprefer(in,'octave',1,octave,0,5);
  in.print  =     xprefer(in,'print',1,0,0,0);
  in.name =       xprefer(in,'name',0,'',0,0);
  in.ftransforms= xcprefer(in,'ftransforms',dgraphs,{zeros(1,maxd)});
  in.gtransforms= xcprefer(in,'gtransforms',dgraphs,in.ftransforms);
  in.minbar =     xcprefer(in,'minbar',dgraphs,{0.01});
  in.scale =      xcprefer(in,'scale',dgraphs,{1.00}); 
  in.cutoff =     xprefer(in,'cutoff',1,-1.e100,0,0);
  in.cutoffs =    xcprefer(in,'cutoffs',dgraphs,{in.cutoff});
  in.graphcutoffs=xcprefer(in,'graphcutoffs',dgraphs,{in.cutoff});
  in.linestyle =  xcprefer(in,'linestyle',dgraphs,{linestyle});
  in.linewidth =  xcprefer(in,'linewidth',dgraphs,{0.5});
  in.legends =    xcprefer(in,'legends',dgraphs,{''});
  in.limits   =   xcprefer(in,'limits',dgraphs,{{zeros(1,2)}});
  in.esample =    xcprefer(in,'esample',dgraphs,{1});
  in.font =       xcprefer(in,'font',dgraphs,{18});
  in.headers =    xcprefer(in,'headers',dgraphs,{''});
  in.logs    =    xcprefer(in,'logs',dgraphs,{{0}});
  in.images =     xcprefer(in,'images',dgraphs,{0});  
  in.imagetype =  xcprefer(in,'imagetype',dgraphs,{1});
  in.slices =     xcprefer(in,'slices',dgraphs,{0});
  in.transverse = xcprefer(in,'transverse',dgraphs,in.slices);
  in.bounds     = xcprefer(in,'bounds',dgraphs,bounds); 
  in.pdimension = xcprefer(in,'pdimension',dgraphs,{3});
  in.diffplot  =  xcprefer(in,'diffplot',dgraphs,{0});
  in.compares =   xcprefer(in,'compares',dgraphs,{0});
  in.gsqplot  =   xcprefer(in,'gsqplot',dgraphs,{0});
  in.chisqplot =  xcprefer(in,'chisqplot',dgraphs,{0});
  in.parametric = xcprefer(in,'parametric',dgraphs,{zeros(1,2)});
  in.compare =    xcprefer(in,'compare',dgraphs,{''});
  in.olabels =    xcprefer(in,'olabels',dgraphs,{' '});
  in.mincount  =  xprefer(in,'mincount',1,10,0,0);
  in.errors =     xprefer(in,'errors',1,0,0,0);  %%Number of error fields
  in.savefig  =   xprefer(in,'savefig',1,0,0,1); %%Saved figures
  in.saveeps  =   xprefer(in,'saveeps',1,0,0,1); %%Saved figures
  if ~isfield(in,'gfunction') && dgraphs > 0     %% If no label data
      for n = in.graphs
        in.gfunction{n} = @(d,~) d;
      end
  end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE GRAPH INDEX
%
  for n = in.graphs                            %% Loop over graphs
    if sequence == 1 
        stn = sprintf(' #%d ',n);
    else
        stn = sprintf(' %d#%d ',s,n);
    end
    if isequal(in.headers{n},' ') || isequal(in.name,' ')
      in.headers{n} = ' ';
    else
      in.headers{n} = [in.name,stn,in.headers{n}];%% Add name to header
    end
    in.xfunctions{n}{maxd}=[];
    in.glabels{n}{maxd}=[];
    in.xk{n}{maxd} = [];
    in.logs{n}{maxd} = [];
    if isfield(in,'xlabels')
        xlabels = in.xlabels;
    end
    if isfield(in,'klabels')
        klabels = in.klabels;
    end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE DIMENSION INDEX
%                          
    for nd = 1:maxd                              %% Loop over dimension  
      if isempty(in.glabels{n}{nd}) && nd<=length(xlabels)
           if nd<=in.dimension && in.gtransforms{n}(nd)
              in.glabels{n}{nd}=klabels{nd};
           else
              in.glabels{n}{nd}=xlabels{nd};
           end
      end                                        %% End if undefined
      if isempty(in.xfunctions{n}{nd})
            in.xfunctions{n}{nd} = @(x,~) x;     %% Return default
      end
      if isempty(in.logs{n}{nd})
            in.logs{n}{nd} = 0;                  %% Return default
      end
    end                                          %% End loop over dimension
  end                                            %% End loop over graphs
  input{s} = in;                                 %% return cells for output
end                                              %% end sequence loop
end                                              %% end function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XGPREFERENCES