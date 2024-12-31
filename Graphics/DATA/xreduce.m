function [d,param,ax,nx,x,xlab,head] =  xreduce(n,d,param,g) 
%   [d,param,ax]  = XREDUCE(n,d,param,g) reduces multidimensional data 
%   for a single xGRAPH plot function to AT MOST three grid dimensions.
%   Input:  graph number 'n', data array 'd', param, parameters 'g'.
%   If g.pdimension{n} < 3 , this defines the maximum number of axes
%   First axis index is time, then x (or kx), y (or ky), etc.
%   Data indices are: line number, axis numbers, error index 
%   Output: reduced data array            'd'
%           reduced parametric array      'param'
%           plotted axis list             'ax'
%           reduced data dimension vector 'nx',
%           reduced coordinate cells      'x', 
%           reduced axis labels           'xlab'
%           reduced observable label      'olab'
%           new header label              'xhead'.
%   Called by: xmultigraph
%   Needs:     xgcompress
%   Licensed by Peter D. Drummond, (2021) - see License.txt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  COMPRESS INPUT DATA 
%    
axes     = g.axes{n};                            %%index vectors of axes
szdi     = size(d);                              %%input data dimensions
leni     = length(szdi);                         %%input data length 
nx       = szdi(1:2);                            %%initialise nx
axsize  = length(szdi)-2;                        %%total axes available
axplots = 0;                                     %%set plot dimensions
ax      = [0,0];
x       = cell(szdi);
xlab    = cell(szdi);
gl      = g.glabels{n};
maxd = min(3,g.pdimension{n});                   %%max dimension plotted
for m = 1:axsize                                 %%loop over dimension
    if szdi(m+1) > 1                             %%if data for plotting
        if axplots >= maxd                       %%if too many axes 
          axes{m} = -1;                          %%set axis point to -1
        else
          if length(axes{m}) > 1                 %if axes points
            axplots = axplots+1;                 %%update axis number
            ax(axplots) = m;                     %%update plot axis list
            nx(axplots+1) = length(axes{m});
            x{axplots} = g.xk{n}{m}(axes{m});
            x{axplots} = g.xfunctions{n}{m}(x{axplots},g);%%update coords
            xlab{axplots} = gl{m};
          end
        end
    end
end
nx(axplots+2) = szdi(leni);
ax(axplots+1) = leni-1; 
g.axes{n} =  axes;
[d,axes]  =  xgcompress(n,d,g);
d   = reshape(d,nx);
if ~isequal(param,0)
    param =  xgcompress(n,param,g);
    param = reshape(param,nx);                   %%reshape parameters
end
head = ' ';
for m = 1:axsize                                 %%loop over dimension
  if isscalar(axes{m})  
    if szdi(m+1) > 1                             %%if data for plotting
        xhead = sprintf('%s=%.3g ',gl{m},g.xk{n}{m}(axes{m}));
        head = [head,xhead];
        xpr(1,g,'Lattice reduced: %s\n',xhead);      
    else
        xpr(1,g,'Lattice reduced: no data for %s axis\n',gl{m});    
    end
  end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XREDUCE FUNCTION