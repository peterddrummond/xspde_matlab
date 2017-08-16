function ec =  xgraph(cdata,input) 
%   ec = XGRAPH(cdata,input) graphs multiple multidimensional data files.
%   Input: data cells 'cdata', input cells 'input'.
%   Output: graphs and sum of comparison differences, `ec'.
%   If no numeric 'cdata' present, reads cdata from a file named 'cdata'.
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

tic();                                         %%set timer
ec = 0;                                        %%Initial comparison errors

%  Get data from file if a filename is present

fileinput = {};                                %%initialize the file input
if ischar(cdata)                               %%If data is character
    [fileinput,cdata] = xread(cdata);          %%Input data from a file 
end                                            %%End if data is character
if nargin == 1                                 %%If no new input data                
    input = fileinput;                         %%get input from a file
end                                            %%End if no new input data 
input = xmakecell(input);                      %%change input to cell
fileinput = xmakecell(fileinput);              %%change file input to cell
input = xgpreferences(input,fileinput);        %%get 'input' defaults
sequence = length(cdata);                      %%get sequence length

%  Loop over a sequence of data inputs

for s = 1:sequence                             %%Loop over sequence
  in = input{s};                               %%inputs for sequence s
  fprintf ('%s sequence %d: %s\n',in.gversion,s,in.name);%%version name
  if in.print >1                               %%if print switch is verbose
    fprintf('xGRAPH data\n');                  %%display lattice data
    display(in);                               %%display lattice data
  end                                          %%end if print switch
  
%  Loop over a set of graph functions
  
  for n =1:in.graphs                           %%Loop over graphs
     if in.pdimension{n} > 0 && in.dimension > 0
       if in.print >2                          %%if graphics print switch 
         fprintf ('Plot %d: %s\n',n,in.gname{n});%%plot name
       end                                    %%end if print switch
       data = cdata{s}{n};
       in.gpoints{n} = size(data);             %%Get n-th input data size 
       data(:,2,:) =  data(:,2,:)+data(:,1,:);
       data(:,3,:) =  data(:,3,:)+data(:,1,:); %%change errors to lines   
       data = real(in.gfunction{n}(data,in));  %%get real graph data
       if  isempty(data)                       %%check if any graph data
             error('xGRAPH error: no data returned from function{%d}\n',n); 
       end                                     %%end check if graph data
       data(:,2,:) =  abs(data(:,2,:)-data(:,1,:));
       data(:,3,:) =  abs(data(:,3,:)-data(:,1,:));%%change lines to errors   
       if  isempty(data)                       %%check if any graph data
           error('No data at function{%d}\n',n); 
       end                                     %%end check if graph data
         
% call multidimensional plot function
 
       err= xmultiplot(n,data,in);             %%single function plot
       if err > 0.0
           fprintf('Max comparison difference in plot %d = %e\n',n,err);
       end
       ec = ec + err;
     end                                       %%end if graph wanted
  end                                          %%end loop over graphs
end                                            %%end sequence
timer = toc();
fprintf('xGRAPH sequence completed, time = %f \n',timer); %%time taken
ec = [ec,timer];
end                                            %%end graphics function