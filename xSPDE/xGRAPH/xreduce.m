function [datan,param,np,x,xlab,olab] =  xreduce(n,datan,param,g) 
%   err = xreduce(n,datan,g) reduces multidimensional data files
%   for a single xSPDE plot function to AT MOST three grid dimensions.
%   Input:  graph number 'n', n-th data array 'datan',parameters 'g'.
%   The g.axes data is used to define the data reduction
%   If: g.axes{n}{nd}<0, the corresponding axis is not plotted
%   If an axis is not plotted, a single default point is used
%   If: g.axes{n}{nd}>=0, the axis is plotted, up to three as a maximum
%   If g.pdimension{n}<3, this defines the maximum number of axes
%   Output: reduced data array            'datan', 
%           reduced data dimension vector 'np', 
%           reduced coordinate cells      'x', 
%           reduced axis labels           'xlab'
%           reduced observable label      'olab'.
%   Licensed by Peter D. Drummond, (2019) - see License.txt 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET PARAMETERS FOR GRAPH (n)
%
    ax =   zeros(1,4);                         %%list of plotted axes    
    np =   g.gpoints{n};                       %%list of points available  
    glab = g.glabels{n};                       %%original x or k label
    olab = g.olabels{n};                       %%vertical axis label
    if g.problength{n} > 1                     %%test for probabilities
         glab{g.dimension+1} = olab;           %%update horizontal labels
         olab  =  ['P(',olab,')'];             %%update vertical axis label
         g.gtransforms{n}=[g.gtransforms{n},0];%%update transform vector
         g.pdimension{n}=[g.pdimension{n},1];  %%update plot dimensions
         g.axes{n}{g.dimension+1} = 1:g.problength{n};%%update axes
         g.kc{g.dimension+1}=0;                %%update k-coordinate
         g.xc{g.dimension+1}=g.probability{n}; %%update x-coordinate
    end                                        %%end test for probabilities
    dimension = length(np)-2;                  %%set data dimensions
    mx = 1+floor(np/2);                        %%index midpoint: space axis
    mx(2) = np(2);                             %%index default: time axis
    maxd = min(3,g.pdimension{n});             %%max dimension plotted
    xlab = glab;                               %%initialize labels
    x = g.xc;                                  %%initialize coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET PLOT AXES, INDEX RANGES
%  
    grd = 1;                                   %%Initial axis index
    for nd = 1:dimension                       %%Loop over dimension
        datan = reshape(datan,np);             %%reshape data to np
        param = reshape(param,np);             %%reshape parameter to np  
        ind = g.axes{n}{nd};                   %%current axis indices
        if length (ind) == 1                   %%if only one index wanted
            if ind == 0                        %%if default index range
                ind = 1:np(grd+1);             %%set full index range
            elseif ind < 0                     %%if no index range
                ind = mx(grd+1);               %%not plotted, set default
            end                                %%end if full index range
        end                                    %%end if one index wanted
        if   grd > maxd                        %%if too many dimensions
            ind = mx(grd+1);                   %%not plotted, set default
        end                                    %%end if full index range
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  REDUCE DATA FOR PART PLOTS
%  
        length_ind =length (ind);              %%points plotted on axis
        n1 = prod(np(1:grd));                  %%total points below
        n2 = np((grd+1));                      %%data points on axis
        n3 = prod(np((grd+2):end));            %%total points above
        datan= reshape(datan,[n1,n2,n3]);      %%reshape data
        param= reshape(param,[n1,n2,n3]);      %%reshape parameters       
        if length_ind  > 1                     %%if axis is plotted
            ax(grd) = nd;                      %%update plotted axis list
            xlab{grd}=glab{nd};                %%update plotted labels
            if g.gtransforms{n}(nd) 
                x{grd} =  g.kc{nd}(ind);       %%graphics k coords
            else
                x{grd} =  g.xc{nd}(ind);       %%get graphics k coords
            end                                %%end if transform switch
            x{grd}=g.xfunctions{n}{nd}(x{grd},g);%%update coordinates
            np(grd+1) = length_ind;            %%update point list
            grd = grd+1;                       %%update axis array index
        else                                   %%too many graphs
            np = [np(1:grd),np((grd+2):end)];
            mx = [mx(1:grd),mx((grd+2):end)];
            xg = g.xc{nd}(ind);
            fprintf('Graph generated at %s = %12.4e\n',glab{nd},xg); 
        end
        datan = datan(:,ind,:);                %%reduce points as required      
        param =param(:,ind,:);                 %%reduce parameter points
    end                                        %%end loop over dimension
    datan = reshape(datan,np);                 %%reshape data
    param = reshape(param,np);                 %%reshape parameters
end