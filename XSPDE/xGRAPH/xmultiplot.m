function err =  xmultiplot(n,datan,g) 
%   err = xmultiplot(n,data,in) plots multidimensional data files
%   for a single xSPDE plot function.
%   Input:  graph number 'n', n-th data array 'datan',parameters 'g'.
%   Output: graphs and maximum comparison differences, `err'.
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

    err = 0;                                   %%Initial comparison errors 
    %set(groot,'defaultAxesFontSize',g.font{n});%%graphics font-size
    set (0, 'defaultaxesfontsize', g.font{n})
    set (0, 'defaulttextfontsize', g.font{n}) 
    olab = g.olabels{n};                       %%get observable labels
    head = g.headers{n};                       %%set graph name
    [datan,nx,x,xlab] =  xreduce(n,datan,g);   %%reduce data dimension
    grd = length(nx)-2;                        %%adjust graph dimension
    if grd == 0
        fprintf('xGRAPH warning: nothing to plot in function{%d}\n',n);
        return
    end
    mx = 1+floor(nx/2);                        %%adjust midpoint: space 
    
    % plot data vs 3D grid 
    
    if grd > 2                             
        ximageplot(n,datan,nx,x,xlab,olab,g);  %%Image plot
        datan=datan(:,:,:,:,mx(5));            %%Central point plotted
    end                                        %%End check if dimension > 2
    
    % plot data vs 2D grid 
    
    if grd  > 1                                %%Check plot dimension
        xtransverseplot(n,datan,nx,x,xlab,olab,g)%%Transverse plot 
        figure;                                %%Start new graphics figure
        im = reshape(datan(1,1,:,:),nx(3),nx(4));%%Image vs t,x
        surf(x{1},x{2},im');                   %%plot 3d mean vs t,x
        xheader(head,xlab{1},xlab{2},olab);    %%3d plot title
        datan=datan(:,:,:,mx(4));              %%Reduce matrix dimension
    end                                        %%End check if dimension > 1
        
    % plot data vs 1D grid      
    xplot(x{1},datan,g.esample{n},g.minbar{n},g.lines{n});
    xheader(head,xlab{1},olab,' ');            %%xspde 2D plot title 
         
    % plot comparison function data vs 1D grid 
    
      
    if ~isempty(g.compare{n})                  %%If comparison results
        da_x = g.compare{n}(x{1},g);           %%get comparison results
        for j = 1:nx(1)
            da(1,1,:) = da_x(j,:);
            plot(x{1},da_x(j,:),'--');         %%2D plot, compare, dashed
            datan(j,1,:) = datan(j,1,:) - da;  %%Store difference
        end
        err= max(max(abs(datan(:,1,:))));      %%Diff vs compare result (n)
        ylabel = strcat('\Delta', olab);       %%2D error plot label
        xplot(x{1},datan,g.esample{n},g.minbar{n},g.lines{n});
        xheader(head,xlab{1},ylabel,' ');      %%Set 2D plot title
    end                                        %%End if compare results
end                                            %%end xplotm function
