    function input = xgpreferences (input,oldinput)  
%   input  =  XGPREFERENCES(input) sets default values for  graphics inputs.
%   Input:  'input' cell array; previous input 'oldinput'.
%   Output: 'input' cell array with updated and default graphics values. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  COMPARE NEW AND OLD SEQUENCES
%                          
if ~isempty(oldinput)                            %% Check  oldinput ~empty
    oldsequence = length(oldinput);              %%get sequence length
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
input = xpreferences(input,'');                  %%get any default inputs                          
sequence = length(input);                        %%get sequence length
lines = {'-k','--k',':k','-.k','-ok','--ok',':ok','-.ok','-+k','--+k'};
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE WHOLE SEQUENCE
%                          
for s = 1:sequence                               %%loop over sequence 
  in = input{s};                                 %%current inputs 
  nd=in.dimension;                               %%number of dimensions 
  xlabels=cell(1,nd);                            %%cell array for xlabels 
  klabels=cell(1,nd);                            %%cell array for klabels
  axes{1} = num2cell(zeros(1,nd));               %%initialize default axes  
  if in.numberaxis || in.dimension > 4
      for i=1:nd
          xlabels{i} = sprintf('r_%d',i);
          klabels{i} = sprintf('k_%d',i);
       end 
  else
      xlabels={'t','x','y','z'};
      klabels={'\omega','k_x','k_y','k_z'};
  end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET COMMON GRAPHICS DEFAULTS
%                          
  in.gversion =   xprefer(in,'gversion',0,'xGRAPH3.1');
  in.graphs =     xprefer(in,'graphs',1,in.functions);
  in.gtransforms = xcprefer(in,'ftransforms',in.graphs,{zeros(1,nd+1)});
  in.axes =       xcprefer(in,'axes',in.graphs,axes);
  in.minbar =     xcprefer(in,'minbar',in.graphs,{0.01});
  in.lines =      xcprefer(in,'lines',in.graphs,{lines});
  in.esample =    xcprefer(in,'esample',in.graphs,{1});
  in.font =       xcprefer(in,'font',in.graphs,{18});
  in.headers =    xcprefer(in,'headers',in.graphs,{in.name});
  in.images =     xcprefer(in,'images',in.graphs,{0});  
  in.imagetype =  xcprefer(in,'imagetype',in.graphs,{1});
  in.transverse = xcprefer(in,'transverse',in.graphs,{0}); 
  in.pdimension = xcprefer(in,'pdimension',in.graphs,{3});
  in.compare =    xcprefer(in,'compare',in.graphs,{''});
  in.parametric = xcprefer(in,'parametric',in.graphs,{zeros(1,2)});  
  if ~isfield(in,'gfunction') && in.graphs>0     %% If no label data
      in.gfunction{in.graphs} = [];
  end
  npoints=prod(in.points);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE GRAPH INDEX
%                          
  for n = 1:in.graphs                            %% Loop over graphs
    if  isempty(in.gfunction{n})
      in.gfunction{n} = @(d,~) d;                %% Return default average
    end                                          %% End if undefined
    in.xfunctions{n}{in.dimension+2}=[];
    in.glabels{n}{in.dimension+2}=[];
    if isfield(in,'xlabels')
        xlabels = in.xlabels;
    end
    if isfield(in,'klabels')
        klabels = in.klabels;
    end
    data = in.gfunction{n}(ones(1,in.errors,npoints),in);
    xfcheck('gfunction',n,data,[1,in.errors,npoints]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER THE DIMENSION INDEX
%                          
    for nd = 1:in.dimension+1                    %% Loop over dimension
      if isempty(in.glabels{n}{nd}) && nd<=length(xlabels)
           if nd<=in.dimension && in.gtransforms{n}(nd)
              in.glabels{n}{nd}=klabels{nd};
           else
              in.glabels{n}{nd}=xlabels{nd};
           end
      end                                        %% End if undefined  
      if isempty(in.xfunctions{n}{nd})
          in.xfunctions{n}{nd} = @(x,~) x;       %% Return default
      else                                       %% End if undefined
          xc = in.xfunctions{n}{nd}(in.xc{nd},in);
          xfcheck('xfunctions',[n,nd],xc,[1,in.points(nd)]);
      end
    end                                          %% End loop over dimension
  end                                            %% End loop over data 
  if in.print > 2
      display(in,'graph parameters');            %% show graph parameters
  end
  input{s} = in;                                 %% return cells for output
end                                              %% end sequence loop
end                                              %% end function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XGPREFERENCES