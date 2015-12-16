function ec =  xgraph(cdata,input) 
%   ec = XGRAPH(cdata,input) graphs multidimensional data files.
%   Input:   data cells 'cdata', input cells 'input'.
%   Output: graphs and max comparison difference, `ec'.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

ec = 0;                                        %%Initialize max errors
if nargin == 2                                 %%If two inputs present
    if ~iscell(input)                          %%if 'input' not cell
        input = {input};                       %%change data to cell
    end                                        %%end if 'input' not cell
    [cdata,input,raw] =  xread(cdata,input);   %%get data if in a file
else                                           %%If one input present
    [cdata,input,raw] =  xread(cdata,{});      %%get data if in a file
end                                            %%end check if 
if ~iscell(cdata)                              %%check if data is not cell
    cdata = {cdata};                           %%change data to cell
end                                            %%end check if 
input = xgrpreferences(input);                 %%get 'input' defaults
sequence = length(cdata);                      %%get sequence length
for s = 1:sequence                             %%Loop over sequence
  data = cdata{s};                             %%Get data cell
  in = input{s};                               %%inputs for sequence s
  nxplot = in.points(1);                       %%time points
  nx = in.points(1:in.dimension);              %%index ranges
  mx = 1+floor(nx/2);                          %%index midpoint for plots
  for n =1:in.graphs                           %%Loop over observables
    set (groot,'defaultAxesFontSize',in.font{n}); %%graphics font-size
    in.images{n} = min(nxplot,in.images{n});   %%image numbers
    in.transverse{n} = min(nxplot,in.transverse{n});
    in.pdimension{n} = min(in.dimension,in.pdimension{n});%%xplot dimension
    head= ' ';                                 %%make graph header blank
    imhead= ' ';                               %%make image header blank
    lab = in.xlabels;                          %%get x labels
    olab = in.olabels{n};                      %%get o labels
    x = in.xc;                                 %%get x coords
    for d = 1:in.pdimension{n}                 %%Loop over dimension
        if in.transforms{n}(d)                 %%If transform switch
            lab{d} = in.klabels{d};            %%get k label
             x{d} = in.gk{d};                  %%get k coords
        end                                    %%End if transform switch
    end                                        %%End dimension loop
    t=x{1};
    compares=0;                                %%comparisons switch
    if ~isempty(in.compare{n})                 %%If comparisons
        compares=1;                            %%If comparisons required
        da_x = in.compare{n}(t,in);            %%get comparison results
    end                                        %%End if comparisons    
    if in.headers{n}                           %%if header wanted 
        head = strcat(in.name,head);           %%set graph name
    end                                        %%end if headers
    datan = reshape(data(:,:,:,n),[3,nx]);     %%extract n-th data      
    if in.dimension  > 3                       %%Check if dimension > 3
        datan = datan(:,:,:,:,mx(4));          %%Central point xplotted
    end                                        %%End if  dimension > 3
    if in.dimension  > 2                       %%Check data dimension > 2
        if in.pdimension{n} > 2                %%Check THIS xplot  dimension
            ximage_xplot()                     %%Image xplot 
        end                                    %%End if g.pdimension > 2
        datan=datan(:,:,:,mx(3));              %%Central point xplotted
    end                                        %%End check if dimension > 2
    if in.dimension  > 1                       %%Check if dimension > 1
      if in.pdimension{n} > 1                  %%Check xplot dimension
        xtransverse_plot()                     %%Transverse xplot 
        figure;                                %%Start new graphics figure
        im = reshape(datan(1,:,:),nx(1),nx(2));%%Image vs t,x
        surf(t,x{2},im');                       %%xplot 3d mean vs t,x
        xheader(head,lab{1},lab{2},olab);      %%xplot title
      end;                                     %%End if xplot dimension
      datan=datan(:,:,mx(2));                  %%Reduce matrix dimension
    end                                        %%End check if dimension > 1
    da_n =datan(1,:);                          %%Reduce data dimension
    eb_n =datan(2,:);                           %%Reduce error dimension     
    se_n =datan(3,:);                          %%Reduce sampling dimension
    xplot(t,da_n,se_n,eb_n,in.errorchecks,t(1),t(nx(1)),in.minbar{n});
    xheader(head,lab{1},olab,' ');             %%xspde 2D plot title
    if compares                                %%If comparison results
      plot(t,da_x,'--');                       %%2D xplot, compare, dashed
      da_n = da_n -  da_x;                     %%Store difference
      err= max(abs(da_n));                     %%Diff vs compare result (n)
      fprintf('Max difference in %d = %e\n',n,err);%%Print errors
      ec = max(err,ec);                        %%Diff vs compared results
      ylabel = strcat('\Delta', olab);         %%2D error xplot label
      xplot(t,da_n,se_n,eb_n,in.errorchecks,t(1),t(nx(1)),in.minbar{n});
      plot(t,0.0*t,'--');                      %%2D xplot, compare, dashed
      xheader(head,lab{1},ylabel,' ');         %%Set 2D xplot title
    end                                        %%End if compare results
  end                                          %%End loop over graphs
end                                            %%End loop over sequence
fprintf('Maximum comparison differences = %e\n',ec); %%Print max diff

function [] =  ximage_xplot()
%% Makes 3d transverse images at fixed time

for i = 0:in.images{n}-1                       %%Loop over images to xplot
    if in.images{n} == 1                       %%only one image wanted
        np = nxplot;                           %%Print last image 
    else                                       %%several images wanted
        np = (nxplot-1)/(in.images{n}-1);      %%image spacing in time
        np = round(1+i*np);                    %%round to next time step   
    end                                        %%End if images == 1
    im = reshape(datan(1,np,:,:),nx(2),nx(3)); %%Image at time t(np)
    if in.headers{n}                           %%if full xheader wanted
        imhead = strcat(olab,sprintf(', t = %0.3f',t(np)));
    end                                        %%end if xheaders
    figure;                                    %%Start new graphics figure
    switch in.imagetype{n}                     %%Select image type
        case 1                                 %first image type
            surf(x{2},x{3},im');               %%xplot 3d surface vs x,y
        case 2                                 %%Second image type
            contourf(x{2},x{3},im');           %%xplot grayxplot vs x,y 
        case 3                                 %%Third image type
            contour(x{2},x{3},im',10);         %%xplot 3d contour vs x,y
    end                                        %%End select image type
    xheader(imhead,lab{2},lab{3},olab);        %%title
end                                            %%End images loop
end                                            %%End image function

function [] =  xtransverse_plot()
%% Makes 2d transverse xplots at fixed time

for i = 0:in.transverse{n}-1                   %%Loop over transverse xplots
    if in.transverse{n} == 1                   %%only one xplot wanted
        np = nxplot;                           %%Print last xplot
    else                                       %%several xplots wanted
        np = (nxplot-1)/(in.transverse{n}-1);  %%xplot spacing in time
        np = round(1+i*np);                    %%round to next time step   
    end                                        %%End if images == 1
    xp = x{2};                                 %%x-coordinate set to x[2}
    if in.headers{n}                           %%if full xheader wanted
       imhead = strcat(olab,sprintf(', t = %0.3f',t(np)));
    end                                        %%end if xheaders
    da_n =reshape(datan(1,np,:),1,nx(2));      %%Transverse data 
    eb_n =reshape(datan(2,np,:),1,nx(2));      %%Transverse error-bar    
    se_n =reshape(datan(3,np,:),1,nx(2));      %%Transverse sampling-error     
    xplot(xp,da_n,se_n,eb_n,in.errorchecks,xp(1),xp(nx(2)),in.minbar{n});
    xheader(imhead,lab{2},olab,'');            %%title
end                                            %%End transverse xplot loop
end                                            %%end Transverse function  
end                                            %%end graphics function


function [] =  xheader(h1,x1,x2,x3)
%% Sets the graph xheader and axis labels

title(h1);                                     %%Set title
xlabel(x1);                                    %%Set x-axis label
ylabel(x2);                                    %%Set y-axis label
zlabel(x3);                                    %%Set z-axis label
end                                            %%end start xplot

function [] =  xplot(t,da_n,se_n,eb_n,errorchecks,lim1,lim2,minbar)
%% Makes 2d xplots with error-bars 

figure;                                        %%Start new graphics figure
relative_error = max(eb_n)/(max(da_n)-min(da_n));
if (errorchecks > 1) && (relative_error > minbar) %%Are there error bars?
    errorbar(t,da_n+se_n,eb_n,'k');            %%Error-bars + sampling
    hold on;                                   %%Allows combined figures
    errorbar(t,da_n-se_n,eb_n,'k');            %%Error-bars - sampling
else                                           %%No error bars!
    plot(t,da_n+se_n,'k');                     %%2D xplot, black, upper sd
    hold on;                                   %%Allows combined figures
    plot(t,da_n-se_n,'k');                     %%2D xplot, black, lower sd
end                                            %%End if errorchecks
xlim([lim1,lim2]);                             %%Set xplot limits
end                                            %%end start xplot



%%Change history:v0.41 
%% - corrected `hold on' position to ensure figures with four sides
%% - test for compares before calling comparison function
%% - add graph number for reference purposes
%% - Change history: v0.42 
%% - added transverse 2d xplot capability
%% - Change history: v0.5
%% - Added data inputs with six dimensions
%% - Change history: v0.51
%% - Modified lattice structures included, added comparison to main 2D graphs
%% - Change history: v0.6
%% - Sequences included