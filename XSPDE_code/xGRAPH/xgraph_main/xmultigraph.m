function [] =  xmultigraph(n,datan,param,g) 
%   [] = XMULTIGRAPH(n,datan,param,g) plots multidimensional data files
%   The first or time axis is plotted from data in 'g.axis' 
%   Input: number 'n', data 'datan', parametric data 'param', structure 'g'.
%   Output: graphs, with multiple plots giving different dimensional views.
%   Called by: xgraph
%   Needs:  xreduce, ximages, xtransverse, xplot3D, xplot2D, xheader
%   xGRAPH functions are licensed by Peter D. Drummond, (2021) - see License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  REDUCE THE INPUT DATA TO 3D
%                          
set (0, 'defaultaxesfontsize', g.font{n})
set (0, 'defaulttextfontsize', g.font{n})     
head = g.headers{n};                             %%set graph name
[datan,param,ax,nx,x,xlab,xhead] =  xreduce(n,datan,param,g); %%reduce data
grd = length(nx)-2;                              %%graph dimension = grd
tlabel = xlab{1};
tcoord = x{1};
olab = g.olabels{n};
if ~isequal(xhead,' ')
    head = [head,xhead];
end
if grd == 0
        xpr(1,g,'xGRAPH warning: nothing to plot in function{%d}\n',n);
        return
end
mx = 1+floor(nx/2);                              %%adjust midpoint 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT DATA IN 3D
%
if grd > 2                                       %Check if dimension > 2
    ximages(n,datan,nx,x,xlab,olab,ax,g);        %%Make images
    sz = size(datan);
    xhead = sprintf('%s=%.3g ',xlab{3},x{3}(mx(4)));
    head = [head,xhead];
    xpr(1,g,'Lattice reduced: %s\n',xhead); 
    datan=reshape(datan(:,:,:,mx(4),:),[sz(1:3),sz(5)]);%%Central point
    if ~isequal(param,0)
        sz = size(param);
        param=reshape(param(:,:,:,mx(4),:),[sz(1:3),sz(5)]);
    end 
end                                              %%End check if dimension > 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT DATA IN 2D
%        
if grd  > 1                                      %%Check plot dimension >1
    xtransverse(n,datan,param,nx,x,xlab,olab,ax(2:end),g) %%Transverse
    im = reshape(datan(1,:,:,1),nx(2),nx(3));    %%Image vs t,x  
    if g.print > 0
        xpr(1,g,'3D plot in graph %d\n',n);
    end
    xplot3D(x{1},x{2},im',xlab{1},xlab{2},olab,head,n,ax,g);
    sz = size(datan);
    xhead = sprintf('%s=%.3g ',xlab{2},x{2}(mx(3)));
    head = [head,xhead];
    xpr(1,g,'Lattice reduced: %s\n',xhead); 
    datan=reshape(datan(:,:,mx(3),:),[sz(1:2),sz(4)]);%%Reduce
    if ~isequal(param,0)
        sz = size(param);
        param=reshape(param(:,:,mx(3),:),[sz(1:2),sz(4)]);
    end
end                                              %%End check if dimension > 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT DATA IN 1D
%            
if g.parametric{n}(2) == 1 && ~isequal(param,0)
    tlabel = g.olabels{g.parametric{n}(1)};
    tcoord = param ;                             %%Axis data 
end
if g.print > 0
    xpr(1,g,'2D plot in graph %d\n',n);
end
xplot2D(tcoord,datan,n,ax,g);                    %%xgraph 2D plot  
xheader(head,tlabel,olab,' ');                   %%xgraph 2D plot title          
end                                              %%end function
