function ec =  xgraph(input,cdata) 
%   ec = XGRAPH(input,cdata) graphs multidimensional data files.
%   Input:  input cells 'input', data cells 'cdata'.
%   Output: graphs and max comparison difference, `ec'.
%   If no numeric 'cdata' present, reads data from a file
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

ec = 0;                                        %%Initialize max errors
if nargin == 2 && ischar(cdata)                %%If data is character
  [hflag,~,cdata,raw] = xread(cdata);          %%get data from a file
  if hflag < 0                                 %%check if file errors
    return                                     %%return
  end                                          %%end check if error
end                                            %%End if data is character
if nargin == 1 || ischar(input)                %%If no input data
  [hflag,input,cdata,raw] = xread(input);      %%get data from a file
  if hflag < 0                                 %%check if file errors
    return                                     %%return
  end                                          %%end check if error
end                                            %%end check if nargin
input = xmakecell(input);                      %%change data to cell
cdata = xmakecell(cdata);                      %%change data to cell
input = xgrpreferences(input);                 %%get 'input' defaults
sequence = length(cdata);                      %%get sequence length

for s = 1:sequence                             %%Loop over sequence  
  in = input{s};                               %%inputs for sequence s
  fprintf ('\nxGRAPH%s: %s, sequence %d\n',in.gversion,in.name,s);%%version
  nxplot = in.points(1);                       %%time points
  nx = in.points(1:in.dimension);              %%index ranges
  mx = 1+floor(nx/2);                          %%index midpoint for plots
  data = reshape(cdata{s},[3,prod(nx),in.graphs]);%%Get data cell
  for n =1:in.graphs                           %%Loop over observables
    set (groot,'defaultAxesFontSize',in.font{n}); %%graphics font-size
    in.images{n} = min(nxplot,in.images{n});   %%image numbers
    in.transverse{n} = min(nxplot,in.transverse{n});
    in.pdimension{n} = min(in.dimension,in.pdimension{n});%%plot dimension
    head= ' ';                                 %%make graph header blank
    imhead= ' ';                               %%make image header blank
    lab = in.xlabels;                          %%get x labels
    olab = in.olabels{n};                      %%get o labels
    x = in.xc;                                 %%get x coords
    if in.headers{n}                           %%if header wanted 
        head = strcat(in.name,head);           %%set graph name
    end                                        %%end if headers
    datan = reshape(data(:,:,n),[3,nx]);       %%extract n-th data   
    for d = 1:in.dimension                     %%Loop over dimension
        if in.transforms{n}(d)                 %%If transform switch
            lab{d} = in.klabels{d};            %%get k label
             x{d} = in.gk{d};                  %%get k coords
             if d > 1                          %%if space transform
                 datan = fftshift(datan,d+1);  %%shift k coords
             end                               %%end if space transform           
        end                                    %%end if transform switch
    end                                        %%end dimension loop
    t=x{1};                                    %%initialize time/frequency
    compares=0;                                %%comparisons switch
    if ~isempty(in.compare{n})                 %%If comparisons
        compares=1;                            %%If comparisons required
        da_x = in.compare{n}(t,in);            %%get comparison results
    end                                        %%End if comparisons    
    if in.dimension  > 3                       %%Check if dimension > 3
        datan = datan(:,:,:,:,mx(4));          %%Central point plotted
    end                                        %%End if  dimension > 3
    if in.dimension  > 2                       %%Check data dimension > 2
        if in.pdimension{n} > 2                %%Check THIS plot  dimension
            ximage_plot()                      %%Image plot 
        end                                    %%End if g.pdimension > 2
        datan=datan(:,:,:,mx(3));              %%Central point plotted
    end                                        %%End check if dimension > 2
    if in.dimension  > 1                       %%Check if dimension > 1
      if in.pdimension{n} > 1                  %%Check plot dimension
        xtransverse_plot()                     %%Transverse plot 
        figure;                                %%Start new graphics figure
        im = reshape(datan(1,:,:),nx(1),nx(2));%%Image vs t,x
        surf(t,x{2},im');                      %%plot 3d mean vs t,x
        xheader(head,lab{1},lab{2},olab);      %%3d plot title
      end;                                     %%End if plot dimension
      datan=datan(:,:,mx(2));                  %%Reduce matrix dimension
    end                                        %%End check if dimension > 1
    da_n =datan(1,:);                          %%Reduce data dimension
    eb_n =datan(2,:);                           %%Reduce error dimension     
    se_n =datan(3,:);                          %%Reduce sampling dimension
    xplot(t,da_n,eb_n,se_n,in.errorchecks,t(1),t(nx(1)),in.minbar{n});
    xheader(head,lab{1},olab,' ');             %%xspde 2D plot title
    if compares                                %%If comparison results
      plot(t,da_x,'--');                       %%2D plot, compare, dashed
      da_n = da_n -  da_x;                     %%Store difference
      err= max(abs(da_n));                     %%Diff vs compare result (n)
      fprintf('Max difference in %d = %e\n',n,err);%%Print errors
      ec = max(err,ec);                        %%Diff vs compared results
      ylabel = strcat('\Delta', olab);         %%2D error plot label
      xplot(t,da_n,eb_n,se_n,in.errorchecks,t(1),t(nx(1)),in.minbar{n});
      plot(t,0.0*t,'--');                      %%2D plot, compare, dashed
      xheader(head,lab{1},ylabel,' ');         %%Set 2D plot title
    end                                        %%End if compare results
  end                                          %%End loop over graphs
end                                            %%End loop over sequence
fprintf('Maximum comparison differences = %e\n',ec); %%Print max diff

function [] =  ximage_plot()
%   XIMAGE_PLOT() makes 3d transverse images at fixed time
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

for i = 0:in.images{n}-1                       %%Loop over images to plot
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
            surf(x{2},x{3},im');               %%plot 3d surface vs x,y
        case 2                                 %%Second image type
            contourf(x{2},x{3},im');           %%plot grayplot vs x,y 
        case 3                                 %%Third image type
            contour(x{2},x{3},im',10);         %%plot 3d contour vs x,y
    end                                        %%End select image type
    xheader(imhead,lab{2},lab{3},olab);        %%title
end                                            %%End images loop
end                                            %%End image function

function [] =  xtransverse_plot()
%   XTRANSVERSE_PLOT() makes 2d transverse images at fixed time
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

for i = 0:in.transverse{n}-1                   %%Loop over transverse plots
    if in.transverse{n} == 1                   %%only one plot wanted
        np = nxplot;                           %%Print last plot
    else                                       %%several plots wanted
        np = (nxplot-1)/(in.transverse{n}-1);  %%plot spacing in time
        np = round(1+i*np);                    %%round to next time step   
    end                                        %%End if images == 1
    xp = x{2};                                 %%x-coordinate set to x[2}
    if in.headers{n}                           %%if full header wanted
       imhead = strcat(olab,sprintf(', t = %0.3f',t(np)));
    end                                        %%end if headers
    da_n =reshape(datan(1,np,:),1,nx(2));      %%Transverse data 
    eb_n =reshape(datan(2,np,:),1,nx(2));      %%Transverse error-bar    
    se_n =reshape(datan(3,np,:),1,nx(2));      %%Transverse sampling-error     
    xplot(xp,da_n,eb_n,se_n,in.errorchecks,xp(1),xp(nx(2)),in.minbar{n});
    xheader(imhead,lab{2},olab,'');            %%title
end                                            %%End transverse plot loop
end                                            %%end Transverse function  
end                                            %%end graphics function

%Version 1.03   data cells are arrays with full space dimensions


function [] =  xheader(h1,x1,x2,x3)
%   XHEADER(h1,x1,x2,x3) makes graph headers and labels axes
%   Here h1 is a header strinng, and x1,x2,x3 are axis label strings
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

title(h1);                                     %%Set title
xlabel(x1);                                    %%Set x-axis label
ylabel(x2);                                    %%Set y-axis label
zlabel(x3);                                    %%Set z-axis label
end                                            %%end xheader

function [] =  xplot(t,da,eb,se,checks,lim1,lim2,minbar)
%   XPLOT(t,da,eb,se,checks,lim1,lim2,minbar)
%   Makes 2D plots with error bars and sampling errors
%   Here 'da' is data, 'se' is sampling error, 'eb' is an error bar
%   Next, 'checks' indicates if error bars should be printed
%   Limits lim1, lim2 are the upper and lower x-axis limits
%   Error bars are only printed if relative size exceeds 'minbar'
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

figure;                                        %%Start new graphics figure
relative_error = max(eb)/(max(da)-min(da));
if (checks > 1) && (relative_error > minbar)   %%Are there error bars?
    errorbar(t,da+se,eb,'k');                  %%Error-bars + sampling
    hold on;                                   %%Allows combined figures
    errorbar(t,da-se,eb,'k');                  %%Error-bars - sampling
else                                           %%No error bars!
    plot(t,da+se,'k');                         %%2D plot, black, upper sd
    hold on;                                   %%Allows combined figures
    plot(t,da-se,'k');                         %%2D plot, black, lower sd
end                                            %%End if errorchecks
xlim([lim1,lim2]);                             %%Set plot limits
end                                            %%end start xplot


