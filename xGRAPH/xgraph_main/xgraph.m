function  xgraph(data,input) 
%   xGRAPH(cdata,input) graphs multidimensional data files.
%   Input: data cells 'data', input parameter cells 'input'.
%   Output: graphs.
%   If no numeric 'data' present, reads data from a file named 'data'.
%   First data dimension is a line index, last dimension indexes error-bars
%   Needs:     xread, xmakecell, xgpreferences, xmultigraph
%   xGRAPH functions are licensed by Peter D. Drummond (2021) - see License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET THE INPUT GRAPH DATA
%
close all;                                       %% Close graphics displays                          
tic();                                           %%set timer
fileinput = {};                                  %%initialize file input
close
if ischar(data)                                  %%If data is character
    [fileinput,data] = xread(data);              %%Input data from a file
    if nargin == 1                               %%If no new input data                
        input = fileinput;                       %%get input from a file
    end                                          %%End if no new input data
else                                             %%End if data is character
    if nargin == 1                               %%If no new input data
        input.name = '';
    end
end
input = xmakecell(input);                        %%change input to cell
fileinput = xmakecell(fileinput);                %%make file input to cell
data = xmakecell(data);                          %%change data to cell
if ~iscell(data{1})                              %%check if data is nested
    data = {data};                               %%create nested cell
end                                              %%end check if nested
sequence = length(data);                         %%get sequence length
for s = 1:sequence                               %%loop over sequence 
    input{s}.dgraphs = length(data{s});          %%get graph numbers
    maxd = 2;                                    %%set data dimension to 2
    sz = cell(1,input{s}.dgraphs);
    for n = 1:input{s}.dgraphs                   %%loop over graphs
        sz{n} = size(data{s}{n});                %%get data size
        maxd = max(maxd,length(sz{n}));          %%get max dimension
    end                                          %%end loop over graphs 
    input{s}.maxd  = maxd;                       %%store max dimension
    input{s}.gpoints = sz;                       %%store data size
end                                              %%end loop over sequence 
input = xgpreferences(input,fileinput);          %%get 'input' defaults
xpr(1,input{1},'\nxGRAPH, v4.0\n');              %%version name
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER INPUT SEQUENCE
%
gtot = 0;
for s = 1:sequence                               %%Loop over sequence
  in = input{s};                                 %%inputs for sequence s
  xpr(1,in,'\nDataset %d: %s\n',s,in.name);      %%version name
  if in.verbose > 1                              %%if print switch verbose
    xpr(2,in,'\nxGRAPH parameters\n');           %%display lattice data
    display(in);                                 %%display lattice data
  end                                            %%end if print switch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOP OVER DATA GRAPHS
%% 
  in.limits{max(in.graphs)+1} = [];              %%Initialize graph limits
  for n = in.graphs                              %%Loop over graphs
    d = data{s}{n};                              %%Get the graph data
    if in.pdimension{n} > 0 && ~isempty(d)       %%Check if graph needed
      xpr(1,in,'\nFunction %d: %s\n',n,in.olabels{n}); %%plot name
      gp = in.gpoints{n};                        %%store input data size
      if in.errors == 0  &&   gp(end) > 1        %%check default case
          gp = [gp,1];                           %#ok<AGROW> %%add index
      end                                        %%end check no errors
      errors = max(1,in.errors);                 %%store corrected errors
      len = length(gp);                          %%store total dimensions
      in.xk{n}{len} = [];                        %%initialize coordinates
      in.axes{n}{len} = [];                      %%initialize axes
      for i = 1:len - 2                          %%loop over data dimension
         if isempty(in.axes{n}{i}) || in.axes{n}{i}(1) == 0 %%if no axes
           in.axes{n}{i} = 1:gp(i+1);            %%use every point
         end                                     %%end check if no axes
         if isempty(in.xk{n}{i})                 %%check if no coordinates 
           in.xk{n}{i} = 1:gp(1+i);              %%use data indices
         end                                     %%end check if coordinates 
      end                                        %%end loop over data
      gpe = gp(end);                             %%last input data size
      compress = prod(gp(1:end-1));              %%compressed data size
      d =  reshape(d,[compress,gpe]);            %%reshape data to matrix
      nc = gpe - errors;                         %%get comparison size
      gp(end) = 3;                               %%set last data dimension
      if nc > 0                                  %%if comparisons exist
         d(:,4:3+nc) = d(:,errors+1:gpe);        %%add comparisons to 4-6  
         d(:,4+nc:6) = 0;                        %%reset comparison errors 
         d(:,5:6) =  d(:,4)+d(:,5:6);            %%get upper error bound
         gp(end) = 6;                            %%set last data dimension
      end                                        %%end if comparisons exist 
      if  errors < 3                             %%if data size<3
         d(:,errors+1:3) = 0;                    %%reset data error fields
      end                                        %%end if data size<3
      d(:,2:3) =  d(:,1) + d(:,2:3);             %%get upper error bounds
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET FUNCTIONS OF DATA GRAPHS
%% 
      d = reshape(d,gp);                         %%reshape data to original
      d = real(in.gfunction{n}(d,in));           %%get real graph function
      if  isempty(d)                             %%check if any graph data
             error('\nNo data in function{%d}\n',n); 
      end                                        %%end check if graph data
      gp =  size(d);                             %%get new data size
      gpe = gp(end);                             %%get new end data size
      gp1 =  gp(1);                              %%get first data size
      ndata  = length(gp)-1;                     %%get new data dimension
      logdat = in.logs{n}{ndata};                %%get log switch for data
      comp = prod(gp(2:end-1));                  %%compressed data size
      d =   reshape(d,[gp1*comp,gpe]);           %%compress data to matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET FUNCTIONAL ERROR DATA 
%% 
      d(:,2:3) =  d(:,2:3) - d(:,1);             %%restore data errors
      if gpe > 3                                 %%if comparisons exist
         if in.chisqplot{n}             
             xchisq(n,reshape(d,gp),in);         %%chi-squared plots
         end
         if in.gsqplot{n}           
             xgsq(n,reshape(d,gp),in);           %%g-squared plots
         end
         dl = zeros([gp1*comp,6]);               %%initialize differences
         d(:,5:6) = d(:,5:6) - d(:,4);           %%restore compare errors
         if logdat                               %%if logs + comparisons
             d(:,1) = max(d(:,1),in.graphcutoffs{n}); %%implement data cutoff
             d(:,4) = max(d(:,4),in.graphcutoffs{n}); %%implement comp cutoff
         end                                     %%end if logs
         gp6 = [gp1,comp,6];                     %%get old array sizes         
         dl(:,1) = d(:,1) - d(:,4);              %%get differences
         dl(:,2:3) = d(:,2:3);                   %%get difference errors
         if max(d(:,5).^2+d(:,6).^2) > 0         %%check if compare error
             dl(:,2) = sqrt(sum(d(:,5:6).^2,2)); %%total compare errors
             dl(:,3) = sqrt(sum(d(:,2:3).^2,2)); %%total data errors
         end                                     %%end check compare error
         sd = sqrt(dl(:,3).^2+dl(:,2).^2);       %%sum the variances
         nc = (sd>0);                            %%valid chi-square
         if in.cutoffs{n}
            nc = nc&(d(:,1)>in.cutoffs{n})&(d(:,4)>in.cutoffs{n});
         end 
         dln = nc.*(dl./(sd+1.e-100));           %%get ratios
         if max(d(:,3).^2+d(:,6).^2) > 0         %%check if sampling error
           cn=sum(nc);                           %%Chi-square valid points
           chs=sum(dln(:,1).^2);                 %%Chi-square error
           if in.verbose > 0                       %%if print switch
             xpr(1,in,'Chi-square errors    = %.3g\n',chs);%%chi-squ error
             xpr(1,in,'Chi-square points    = %d\n',cn);%%comparison points
           end                                   %%end if print switch
         end                                     %%end check sampling error
         if in.diffplot{n} == 2                  %%check if normalized diff
           dl=dln;                               %%normalize the diff
         end                                     %%end check if normalized
         if in.diffplot{n} < 3                   %%if diff-switch < 3     
             dl = reshape(dl,gp6);               %%reshape differenecs     
             dl = cat(1,dl(:,:,1:3),dl(:,:,4:6));%%make diffs into lines
             gp = [2*gp1,gp(2:end-1),3];         %%get new array sizes
             dl = reshape(dl,gp);                %%reshape diffs for plots
             d = reshape(d,gp6);                 %%reshape data 
             d = cat(1,d(:,:,1:3),d(:,:,4:6));   %%make data into lines    
         else                                    %%if diff-switch == 3
             dl =  d(:,4:6);                     %%get pure comparisons
             gp = [gp1,gp(2:end-1),3];           %%get new array sizes 
             dl = reshape(dl,gp);                %reshape diffs for plots
             d =   d(:,1:3) ;                    %%make data into lines   
         end                                     %%end if diff-switch
         d = reshape(d,gp);                      %%reshape data for plots      
      else                                       %%else no comparisons
         d =   reshape(d(:,1:3),gp);             %%reshape data for plots
      end                                        %%end if comparisons exist 
      in.gpoints{n} = gp;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  GET PARAMETRIC AXIS DATA
%% 
      if in.parametric{n}(1) > 0                 %%if parametric plot
         param = data{s}{in.parametric{n}};      %%parametric data included
         sizepa = size(param);
         if ~isequal(sizepa(1),gp(1))            %%if line number mismatch
            error('Graph %d, line mismatch: parametric = %d, data = %d'...
                  ,n,sizepa(1),gp(1));           %%print error message
         end                                     %%end if line mismatch
      else
         param = zeros(gp);
      end                                        %%end if parametric plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CALL MULTIPLE PLOT FUNCTION
%% 
      xmultigraph(n,d,param,in);                 %%multidimensional plot
      if gpe > 3                                 %%if compare & difference
        if in.diffplot{n} > 0 
        in.graphcutoffs{n} = -1.e100;            %%remove graph cutoff
        in.logs{n}{ndata} = 0;                   %%remove log switch
        switch in.diffplot{n}        
          case 1
          in.olabels{n}=['\Delta',in.olabels{n}];%%delta label
          case 2
          in.olabels{n}=['\Delta',in.olabels{n},'/\sigma'];
          otherwise
          in.olabels{n}=['Cf ',in.olabels{n}];   %%pure comparison data
        end                                      %%end if not normalized
        xmultigraph(n,dl,param,in);              %%comparison plots
        end
      end                                        %%end if difference plot
    end                                          %%end if graph wanted
    gtot = gtot + 1;
 end                                             %%end loop over graphs
end                                              %%end sequence
xpr(0,in,'\nxGRAPH plotted %d input(s), time = %.3gs\n',gtot,toc());
h=findobj('type','figure');                      %% find figure handles
for k=1:numel(h)                                 %% loop on handles
    fname=sprintf('Fig%d',k);                    %% Set filename
    if in.savefig                                %% If figure save = Y
         savefig(h(k),[fname,'.fig']);           %% Save fig k
    end                                          %% End if figure save = Y
    if in.saveeps                                %% If EPS save = Y
          saveas(h(k),[fname,'.eps'],'epsc');    %% Save EPS k
    end                                          %% End if EPS save = Y
end                                              %% End loop on handles
end                                              %%end graphics function