function [datan,np,x,xlab] =  xreduce(n,datan,g) 
%   err = xreduce(n,datan,g) reduces multidimensional data files
%   for a single xSPDE plot function in at most three grid dimensions.
%   Input:  graph number 'n', n-th data array 'datan',parameters 'g'.
%   Output: reduced data, points, coordinates, labels.
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 


    grd = 0;                                   %%Initial graph number 
    ax = zeros(1,4);                           %%list of plotted axes
    x = g.xc;                                  %%x coords
    np = g.gpoints{n};                         %%index ranges
    dimension = length(np)-2;
    mx = 1+floor(np/2);                        %%index midpoint: space axis
    mx(3) = np(3);                             %%index default: time axis
    xlab = g.glabels{n};                       %%get x or k label

    %Loop over dimension to find data axes to plot
    
    for nd = 1:dimension                       %%Loop over dimension
        datan = reshape(datan,np);             %%reshape data to np
        if g.transforms{n}(nd)                 %%If transform switch
            x{nd} =  fftshift(g.kc{nd});       %%get graphics k coords
            datan = fftshift(datan,nd+2);      %%shift data k coords        
        end                                    %%end if transform switch
        if length (g.axes{n}{nd}) == 1         %%if only one index wanted
            if g.axes{n}{nd} == 0              %%if default index range
                g.axes{n}{nd} = 1:np(nd+2);    %%set full index range
            elseif g.axes{n}{nd} < 0           %%if negative index range
                g.axes{n}{nd} = mx(nd+2);      %%not plotted, set default
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
                g.axes{n}{nd} = mx(nd+2);      %%set default index 
            end                                %%end if less than 3 graphs
        end                                    %%end if axis plotted
        n1 = prod(np(3:nd+1));                 %%points below
        n2 = prod(np(nd+3:dimension+2));       %%points above
        datan= reshape(datan,[np(1:2),n1,np(nd+2),n2]);%%reshape data
        datan = datan(:,:,:,g.axes{n}{nd},:);  %%reduce points
        np(nd+2) = length(g.axes{n}{nd});      %%reset points count
    end                                        %%end loop over dimension
    np = [np(1:2),np(2+ax(1:grd))];             %%adjust points plotted
    datan = reshape(datan,np);                 %%reshape data
end