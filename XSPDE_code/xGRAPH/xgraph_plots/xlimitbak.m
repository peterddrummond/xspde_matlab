function [] =  xlimit(n,ax,g)
%   XLIMIT(n,ax,g)
%   Uses the Matlab xlim, ylim, zlim function to set plot limits
%   g.limit{n}={[tl,tu],[xl,xu],..} is the n-th plot limit array
%   The current axis list for plot n is denoted ax
%   Last axis limit is the observable, with dimension d+1 for d dimensions
%   Sets either two or three limits, depending on length of ax
%   If plot range from limit data is zero or nonexistent, no limits are set
%   Plots have up to three dimensions
%   xGRAPH functions are licensed by Peter D. Drummond, (2020) - see License

limit = g.limits{n};                    %%Get n-th plot limits 
sz = length(limit);                     %%Get length of plot limits
lx = length(ax);                        %%Get length of axis vector
axis = gca;                             %%Get current axis handle
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FIRST COORDINATE PLOT LIMITS       

if sz >= ax(1) && ax(1) > 0             %%If limit exists   
    lim = limit{ax(1)};                 %%Get coordinate limits                             
    if lim(1) < lim(2)                  %%If plot range nonzero   
        axis.XLim=lim;                  %%Set coordinate limits  
    end                                 %%End if plot range nonzero 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SECOND COORDINATE PLOT LIMITS
        
if sz >= ax(2) && ax(1) > 0             %%If limit exists   
    lim = limit{ax(2)};                 %%Get coordinate limits                             
    if lim(1) < lim(2)                  %%If plot range nonzero   
        axis.YLim=lim;                  %%Set coordinate limits  
    end                                 %%End if plot range nonzero 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   THIRD COORDINATE PLOT LIMITS
        
if lx > 2 && sz >= ax(3) && ax(3) > 0   %%If axis and limit exists   
    lim = limit{ax(3)};                 %%Get coordinate limits                             
    if lim(1) < lim(2)                  %%If plot range nonzero   
        axis.ZLim=lim;                  %%Set coordinate limits  
    end                                 %%End if plot range nonzero 
end                                     %%End if axis and limit exists
end                                     %%End function     