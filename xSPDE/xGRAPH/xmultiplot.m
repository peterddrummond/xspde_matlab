function err =  xmultiplot(n,datan,param,g) 
%   err = xmultiplot(n,data,axis,in) plots multidimensional data files
%   for a single xSPDE plot function.
%   The first or time axis is plotted from data in 'axis' 
%   Input:  number 'n', data & aixs array 'datan', 'axis',parameters 'g'.
%   Output: graphs and maximum comparison differences, `err'.
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

    err = 0;                                   %%Initial comparison errors 
    set (0, 'defaultaxesfontsize', g.font{n})
    set (0, 'defaulttextfontsize', g.font{n})     
    head = g.headers{n};                       %%set graph name   
    [datan,param,nx,x,xlab,olab] =  xreduce(n,datan,param,g); %%reduce data
    grd = length(nx)-2;                        %%adjust graph dimension
    tlabel = xlab{1};
    tcoord = x{1};
    if grd == 0
        fprintf('xGRAPH warning: nothing to plot in function{%d}\n',n);
        return
    end
    mx = 1+floor(nx/2);                        %%adjust midpoint: space 
    
    % plot data vs 3D grid 
    
    if grd > 2                             
        ximageplot(n,datan,nx,x,xlab,olab,g);  %%Image plot
        sz = size(datan);
        datan=reshape(datan(:,:,:,mx(4),:),[sz(1:3),sz(5)]);%%Central point
        sz = size(param);
        param=reshape(param(:,:,:,mx(4),:),[sz(1:3),sz(5)]);%%Central point 
    end                                        %%End check if dimension > 2
    
    % plot data vs 2D grid 
    
    if grd  > 1                                %%Check plot dimension
        xtransverseplot(n,datan,param,nx,x,xlab,olab,g)%%Transverse plot 
        figure;                                %%Start new graphics figure
        im = reshape(datan(1,:,:,1),nx(2),nx(3));%%Image vs t,x
        surf(x{1},x{2},im');                   %%plot 3d mean vs t,x
        xheader(head,xlab{1},xlab{2},olab);    %%3d plot title
        sz = size(datan);
        datan=reshape(datan(:,:,mx(3),:),[sz(1:2),sz(4)]);%%Reduce
        sz = size(param);
        param=reshape(param(:,:,mx(3),:),[sz(1:2),sz(4)]);
    end                                        %%End check if dimension > 1
        
    % plot data vs 1D grid
    
    if g.parametric{n}(2) == 1
        tlabel = g.olabels{g.parametric{n}(1)};
        tcoord =param ;%%Axis data 
    end
    xplot(tcoord,datan,g.esample{n},g.minbar{n},g.lines{n});
    xheader(head,tlabel,olab,' ');            %%xspde 2D plot title 
         
    % plot comparison function data vs 1D grid 
    
      
    if ~isempty(g.compare{n})                  %%If comparison results
        da_x = g.compare{n}(x{1},g);           %%get comparison results
        if length(da_x) == 1
            da_x = da_x+0*x{1};
        end
        sc = xfcheck('compare',n,da_x,[0,length(x{1})]);
        datan = datan(1:sc(1),:,:);
        for j = 1:sc(1)
            da(1,:,1) = da_x(j,:);
            plot(tcoord,da_x(j,:),'-.');       %%2D plot, compare, dashed
            datan(j,:,1) = datan(j,:,1) - da;  %%Store difference
        end
        err= max(max(abs(datan(:,:,1))));      %%Diff vs compare result (n)
        ylabel = strcat('\Delta', olab);       %%2D error plot label
        xplot(tcoord,datan,g.esample{n},g.minbar{n},g.lines{n});
        xheader(head,tlabel,ylabel,' ');       %%Set 2D plot title
    end                                        %%End if compare results
end                                            %%end xplotm function
