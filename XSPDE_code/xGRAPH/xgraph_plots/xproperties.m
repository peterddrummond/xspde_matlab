function [] =  xproperties(n,ax,g)
%   XPROPERTIES(n,ax,g)
%   Uses the Matlab graphics functions to set n-th graph properties
%   The current axis list for graph n is denoted ax
%   Here ind = ax(1) is the index of the plotted x-coordinate
%   and ind = ax(end) is the index of the dependent or o-coordinate
%   Note that initially, ind = 1,2,..d+1 for d space-time dimensions
%   First axis limit is the t-axis 
%   g.limits{n}=[tl,tu;xl,xu;,..;ol,ou] is the n-th plot limit matrix
%   Sets either two or three limits, depending on length of ax
%   If plot range from limit data is zero or nonexistent, no limits are set
%   If log switches are set, uses a log scale
%   One, two or all three axes can have log scales
%   Plots have up to three dimensions, chosen out of the data dimensions
%   xGRAPH functions are licensed by Peter D. Drummond, (2020) - see License

limits = g.limits{n};                    %%Get n-th plot limits
logs = g.logs{n};                        %%Get n-th plot log switches
sz1 = length(limits);                    %%Get length of plot limits
lx = length(ax);                         %%Get length of axis vector
axis = gca;                              %%Get current axis handle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FIRST COORDINATE PLOT LIMITS       

if ax(1) > 0                            %%If axis exists      
  if sz1 >= ax(1) && ~isempty(limits)   %%If limit exists   
    lim = limits{ax(1)};                %%Get coordinate limits                             
    if length(lim) == 2 && lim(1) < lim(2)  %%If plot range nonzero   
        axis.XLim=lim;                  %%Set coordinate limits  
    end                                 %%End if plot range nonzero 
  end                                   %%End if limit exists 
  if logs{ax(1)}                        %%If log switch is set  
    axis.XScale='log';                  %%Set log scales  
  end                                   %%End if log switch is set 
end                                     %%End if axis exists
                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SECOND COORDINATE PLOT LIMITS
        
if length(ax)>1 && ax(2) > 0            %%If axis exists      
  if sz1 >= ax(2) && ~isempty(limits)   %If limit exists   
    lim = limits{ax(2)};                %%Get coordinate limits                              
    if length(lim) == 2 && lim(1) < lim(2) %%If plot range nonzero   
        axis.YLim=lim;                  %%Set coordinate limits  
    end                                 %%End if plot range nonzero 
  end                                   %%End if limit exists 
  if logs{ax(2)}                        %%If log switch is set  
    axis.YScale='log';                  %%Set log scales  
  end                                   %%End if log switch is set 
end                                     %%End if axis exists

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   THIRD COORDINATE PLOT LIMITS
        
if length(ax)>2 && ax(3) > 0            %%If axis exists  
  if lx > 2 && sz1 >= ax(3) && ~isempty(limits)    
    lim = limits{ax(3)};                %%Get coordinate limits                              
    if length(lim) == 2 && lim(1) < lim(2) %%If plot range nonzero   
        axis.ZLim=lim;                  %%Set coordinate limits  
    end                                 %%End if plot range nonzero
  end                                   %%End if limit exists 
  if logs{ax(3)}                        %%If log switch is set  
    axis.ZScale='log';                  %%Set log scales  
  end                                   %%End if log switch is set 
end                                     %%End if axis exists
end                                     %%End function     