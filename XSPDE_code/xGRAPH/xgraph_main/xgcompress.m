function [data,axes]  =  xgcompress(n,data,r) 
%   a = XGCOMPRESS(n,data,r) compresses graphics data on the required axes. 
%   Input: index 'n', data array 'data', parameter structure 'r'
%   Output: compressed data with new axes cell arrays.
%   The input axes parameter is used to define the data reduction.
%   The first data dimension is not compressed: it is not an axis dimension
%   If axes{nd} = -2, the last point is used
%   If axes{nd} = -1, the midpoint is used
%   If axes{nd} = 0, all the data is retained
%   If axes{nd} > 0, the specified data point or points are retained
%   Corrected axes and compressed data is returned.
%   Called by: xreduce
%   xGRAPH functions are licensed by Peter D. Drummond, (2021) - see License 
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE PLOTTED AXIS DATA
% 
    np = size(data);                             %%list of points available
    axes = r.axes{n};                            %%list of axis points
    dimension = min(length(axes),length(np)-2);  %%data axis dimensions
    mid = 1+floor(np/2);                         %%midpoint of axis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET PLOT AXES, INDEX RANGES
%  
    for nd = 1:dimension                         %%Loop over dimension
        ind = axes{nd};                          %%axis index list for nd
        if length (ind) == 1                     %%if only one index wanted
            if ind == 0                          %%if default index range
                ind = 1:np(nd+1);                %%set full index range
            elseif ind == -1                     %%if index = -1
                ind = mid(nd+1);                 %%set to midpoint 
            elseif ind == -2                     %%if index = -2
                ind = np(nd+1);                  %%set to last value
            end                                  %%end if full index range
        end                                      %%end if one index
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  REDUCE DATA FOR PARTIAL PLOT
%        
        if np(nd+1) > length(ind)                %%more data than indices
          n1 = prod(np(1:nd));                   %%total points below
          n2 = np((nd+1));                       %%data points on axis
          n3 = prod(np((nd+2):end));             %%total points above
          data= reshape(data,[n1,n2,n3]);        %%reshape data
          np(nd+1) = length(ind);                %%update point list
          data = data(:,ind,:);                  %%reduce data as required
        else                                     %%less data than indices
            ind = ind(1:np(nd+1));               %%reduce index to data
        end                                      %%end if more data
        axes{nd} = ind;                          %%restore axis index list
    end                                          %%end loop over dimension
    data = reshape(data,np);                     %%reshape data    
end                                              %%end xcompress function