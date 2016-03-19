function [datan,nx,x,xlab] =  xreduce(n,datan,g) 
%   err = xreduce(n,datan,g) reduces multidimensional data files
%   for a single xSPDE plot function in at most three grid dimensions.
%   Input:  graph number 'n', n-th data array 'datan',parameters 'g'.
%   Output: reduced data, points, coordinates, labels.
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 


    grd = 0;                                   %%Initial graph number 
    ax = zeros(1,4);                           %%list of plotted axes
    x = g.xc;                                  %%x coords
    nx = g.points;                             %%index ranges
    mx = 1+floor(nx/2);                        %%index midpoint: space axis
    mx(1) = nx(1);                             %%index default: time axis
    xlab = g.glabels{n};                       %%get x or k label
    
    %Loop over dimension to find data axes to plot
    
    for nd = 1:g.dimension                     %%Loop over dimension      
        datan = reshape(datan,[3,nx]);          %%reshape data
        if g.transforms{n}(nd)                 %%If transform switch
            x{nd} =  fftshift(g.kc{nd});       %%get graphics k coords
            datan = fftshift(datan,nd+1);      %%shift data k coords        
        end                                    %%end if transform switch
        if length (g.axes{n}{nd}) == 1         %%if one number in range
            if g.axes{n}{nd} == 0              %%if default index range
                g.axes{n}{nd} = 1:nx(nd);      %%set full index range
            elseif g.axes{n}{nd} < 0           %%if negative index range
                g.axes{n}{nd} = mx(nd);        %%not plotted, set default
            end                                %%end if full index range
        end                                    %%end if one number in range
        if length (g.axes{n}{nd}) > 1          %%if axis has finite range
            if grd < min(3,g.pdimension{n})    %%if less than max dimension
                grd = grd+1;                   %%update number of graphs
                ax(grd) = nd;                  %%update plotted axes
                xlab{grd}=xlab{nd};            %%update plotted labels
                xnd = x{nd}(g.axes{n}{nd});
                x{grd}=g.xfunctions{n}{nd}(xnd,g);%%update coordinates
            else                               %%too many graphs
                g.axes{n}{nd} = mx(nd);        %%set default index 
            end                                %%end if less than 3 graphs
        end                                    %%end if axis plotted
        n1 = prod(nx(1:nd-1));                 %%points below
        n2 = prod(nx(nd+1:g.dimension));       %%points above
        datan= reshape(datan,[3,n1,nx(nd),n2]);%%reshape data
        datan = datan(:,:,g.axes{n}{nd},:);    %%reduce points
        nx(nd) = length(g.axes{n}{nd});        %%reset points count   
    end                                        %%end loop over dimension
    nx = nx(ax(1:grd));                        %%adjust points plotted
    if ~isempty(nx)
        datan = reshape(datan,[3,nx]);         %%reshape data
    end
end