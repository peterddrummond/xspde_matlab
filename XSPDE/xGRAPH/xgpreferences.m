function input = xgpreferences (input,oldinput)  
%   input  =  XGRPREFERENCES(input) sets default values for  graphics inputs.
%   Input:  'input' cell array. previous input 'oldinput'.
%   Output: 'input' cell array with updated and default graphics values. 
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

                          %% Compare new and old input sequence 

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

                          %% Obtain the number of graph functions 
                          
sequence = length(input);                        %%get sequence length
for s = 1:sequence                               %%loop over sequence 
  in = xpreferences(input{s});                   %%get input structure
  graphs=0;                                      %%Initial plot number
  if isfield(in,'function')                      %%If plot functions input
      graphs = length(in.function);              %%Number of plot functions
  end                                            %%End if plot functions 
  if graphs < in.averages                        %%More averages than plots
      graphs = in.averages;                      %%Plots set to averages
      in.function{graphs} =[];                   %%Set last plot function
  end                                            %%End averages vs plots

                       %% Set the graphics default parameters
                       
  if in.ensembles(2)*in.ensembles(3) > 1         %%If sampling errors
      esample = 1;                               %%Set sampling flag to 1
  else                                           %%Else no sampling errors
      esample = 0;                               %%Set sampling flag to 0
  end                                            %%End if sampling errors
  if in.numberaxis || in.dimension > 4
      xlabels={'t','x_1','x_2','x_3','x_4','x_5','x_6','x_7','x_8'};
      klabels={'\omega','k_1','k_2','k_3','k_4','k_5','k_6','k_7','k_8'};
  else
      xlabels={'t','x','y','z'};
      klabels={'\omega','k_x','k_y','k_z'};
  end
  
  in.gversion =   xprefer(in,'gversion',0,'xGRAPH1.1');
  in.graphs =     xprefer(in,'graphs',1,graphs);
  in.transforms = xcprefer(in,'transforms',in.graphs,{zeros(1,in.dimension)});
  in.axes =       xcprefer(in,'axes',in.graphs,{num2cell(zeros(1,in.dimension))});
  in.minbar =     xcprefer(in,'minbar',in.graphs,{0.01});
  in.ebar =       xcprefer(in,'ebar',in.graphs,{in.ebar});
  in.esample =    xcprefer(in,'esample',in.graphs,{esample});
  in.font =       xcprefer(in,'font',in.graphs,{18});
  in.headers =    xcprefer(in,'headers',in.graphs,{in.name});
  in.images =     xcprefer(in,'images',in.graphs,{0});  
  in.imagetype =  xcprefer(in,'imagetype',in.graphs,{1});
  in.transverse = xcprefer(in,'transverse',in.graphs,{0}); 
  in.pdimension = xcprefer(in,'pdimension',in.graphs,{3});
  in.compare =    xcprefer(in,'compare',in.graphs,{''});
  
  for n = 1:graphs                               %% Loop over graphs
    if  isempty(in.function{n})
      if  n<=in.averages
          in.function{n} = @(d,~) d{n}(:,:,:);   %% Return default average
      else
          error ('xGRAPH error: no function, sequence %d, graph %d\n',s,n);
      end
    end                                          %% End if undefined
    in.xfunctions{n}{in.dimension+1}=[];
    in.glabels{n}{in.dimension+1}=[];
    if isfield(in,'xlabels')
        xlabels = in.xlabels;
    end
    if isfield(in,'klabels')
        klabels = in.klabels;
    end
    for nd = 1:in.dimension                      %% Loop over dimension
       if isempty(in.glabels{n}{nd})
           if in.transforms{n}(nd)
              in.glabels{n}{nd}=klabels{nd};
           else
              in.glabels{n}{nd}=xlabels{nd};
           end
      end                                        %% End if undefined  
      if isempty(in.xfunctions{n}{nd})
          in.xfunctions{n}{nd} = @(x,~) x;       %% Return default
      end                                        %% End if undefined
    end                                          %% End loop over dimension
  end                                            %% End loop over data 
  
  if in.print > 2
      display(in,'graph parameters');            %% show graph parameters
  end
  input{s} = in;                                 %%return cells
end                                              %% sequence loop
end                                              %% end function