function err =  xmultiplot(n,datan,g) 
%   err = xmultiplot(n,data,in) plots multidimensional data files
%   for a single xSPDE plot function.
%   Input:  graph number 'n', n-th data array 'datan',parameters 'g'.
%   Output: graphs and maximum comparison differences, `err'.
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

    err = 0;                                   %%Initial comparison errors 
    set(groot,'defaultAxesFontSize',g.font{n});%%graphics font-size   
    olab = g.olabels{n};                       %%get observable labels
    head = g.headers{n};                       %%set graph name
    [datan,nx,x,xlab] =  xreduce(n,datan,g);  %%reduce data dimension
    grd = length(nx);                          %%adjust graph dimension
    if grd == 0
        fprintf('xGRAPH warning: nothing to plot in function{%d}\n',n);
        return
    end
    mx = 1+floor(nx/2);                        %%adjust midpoint: space 
    mx(1) = nx(1);                             %%adjust default: time 
    
    % plot data vs 3D grid 
    
    if grd > 2                             
        ximageplot(n,datan,nx,x,xlab,olab,g); %%Image plot
        datan=datan(:,:,:,mx(3));              %%Central point plotted
    end                                        %%End check if dimension > 2
    
    % plot data vs 2D grid 
    
    if grd  > 1                                %%Check plot dimension
        xtransverseplot(n,datan,nx,x,xlab,olab,g)%%Transverse plot 
        figure;                                %%Start new graphics figure
        im = reshape(datan(1,:,:),nx(1),nx(2));%%Image vs t,x
        surf(x{1},x{2},im');                   %%plot 3d mean vs t,x
        xheader(head,xlab{1},xlab{2},olab);   %%3d plot title
        datan=datan(:,:,mx(2));                %%Reduce matrix dimension
    end                                        %%End check if dimension > 1
        
    % plot data vs 1D grid 
      
    da_n =datan(1,:);                          %%Reduce data dimension
    eb_n =datan(2,:);                          %%Reduce error dimension
    if g.esample{n}                            %%If sampling error required
        se_n =datan(3,:);                      %%Get sampling error  
    else                                       %%No sampling error exists
        se_n = zeros(1,nx(1));                 %Set sampling error to zero 
    end                                        %%End if sampling error 
    x2plot(x{1},da_n,eb_n,se_n,g.ebar{n},g.minbar{n});
    xheader(head,xlab{1},olab,' ');            %%xspde 2D plot title 
         
    % plot comparison function data vs 1D grid 
      
    if ~isempty(g.compare{n})                  %%If comparison results
        da_x = g.compare{n}(x{1},g);           %%get comparison results
        plot(x{1},da_x,'--');                  %%2D plot, compare, dashed
        da_n = da_n -  da_x;                   %%Store difference
        err= max(abs(da_n));                   %%Diff vs compare result (n)
        ylabel = strcat('\Delta', olab);       %%2D error plot label
        x2plot(x{1},da_n,eb_n,se_n,g.ebar{n},g.minbar{n});
        plot(x{1},0.0*x{1},'--');              %%2D plot, compare, dashed
        xheader(head,xlab{1},ylabel,' ');      %%Set 2D plot title
    end                                        %%End if compare results
end                                            %%end xplotm function
