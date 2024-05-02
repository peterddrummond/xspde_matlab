function [data,x,axes]  =  xcompress(n,data,p) 
%   a = XCOMPRESS(n,data,p) compresses the data on the required axes. 
%   Input: index 'n', data array 'data', parameter structure 'p'
%   Output: compressed data with new coordinate and axes cell arrays.
%   The axes parameter is used to define the data reduction, as follows:
%   Axes: 'p.axes'; x,k-coordinates: 'p.xc','p.kc', switch: 'p.ftransforms' 
%   The first data dimension is not compressed: it is not an axis dimension
%   If axes{nd}=-2, the last point is used
%   If axes{nd}=-1, the midpoint is used
%   If axes{nd}=0, all the data is retained
%   If axes{nd}>0, the specified data point or points are retained
%   Called by xsim and used for comparison and chi-square error estimates
%   xSPDE functions are licensed by Peter D. Drummond, (2021) - see License 
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE PLOTTED AXIS DATA
% 
np = size(data);                                 %%list of points available
axes = p.axes{n};                                %%list of axis points
mid = 1+floor(np/2);                             %%midpoint of axis
x =   p.xc;                                      %%initialize coordinates
m1 = length(x);                                  %%length of coord. list
if isfield(p,'oc')                               %%check if bin coordinates
    if n <= length(p.oc)                         %%check if n-th bin
        oc = p.oc{n};                            %%n-th bin coordinates
        for m = 1:length(oc)                     %%loop over bin coords 
            x{m+m1} = oc{m};                     %%store bin coordinates
        end                                      %%end loop on bins
    end                                          %%end check if n-th bin
end                                              %%end check if bins
d = length(p.gpoints{n})-2;                      %%data axis dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET PLOT AXES, INDEX RANGES
%  
for nd = 1:d                                     %%Loop over dimension   
    ind = axes{nd};                              %%axis index list  
    xnd = x{nd};                                 %%axis nd coordinates
    if nd <= p.dimensions && p.ftransforms{n}(nd)%%if transforms are set
        xnd =  p.kc{nd};                         %%graphics k coords     
    end                                          %%end if space dimension
    if isscalar(ind)                             %%if only one index wanted
        if ind == 0                              %%if default index range
            ind = 1:np(nd+1);                    %%set full index range
        elseif ind == -1                         %%if index = -1
            ind = mid(nd+1);                     %%set to midpoint 
        elseif ind == -2                         %%if index = -2
            ind = np(nd+1);                      %%set to last value
        end                                      %%end if full index range
    end                                          %%end if one index
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  REDUCE DATA FOR PARTIAL PLOT
%  
   if np(nd+1) > length(ind)                     %%if more data than needed
        n1 = prod(np(1:nd));                     %%total points below
        n2 = np((nd+1));                         %%data points on axis
        n3 = prod(np((nd+2):end));               %%total points above
        data= reshape(data,[n1,n2,n3]);          %%reshape data
        np(nd+1) = length(ind);                  %%update point list
        data = data(:,ind,:);                    %%reduce data as required
    else                                         %%if less data than needed
        ind = ind(1:np(nd+1));                   %%reduce index to data
    end                                          %%end if more data
    x{nd} = xnd(ind);                            %reset coordinates
    axes{nd} = 1:length(x{nd});                  %reset axes
end                                              %%end loop over dimension
data = reshape(data,np);                         %%reshape data    
end                                              %%end xcompress function