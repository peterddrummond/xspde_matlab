function [] =  xplot2D(t,dat,n,ax,g)
%   XPLOTD(t,dat,n,ax,g)
%   Makes 2D plots with error bars and sampling errors against axis 't'
%   Here 'dat' is the n-th graph data, with sampling error and error bars
%   The ax input determines which axes are plotted.
%   Note: 't' is a dummy label for axis-1, it could be anything
%   Input data has three dimensions: lines, points and errors
%   Input t has up to three dimensions, must have at least one
%   The value of 'esample' indicates how sampling error is graphed
%   If g.esample{n} = 0, no sampling error lines are plotted, just the mean
%   If g.esample{n} = m, m*sigma sampling errors are included in errorbars
%   If g.esample{n} = -1, combined sum error bars plotted
%   If g.esample{n} = -2, combined rms error bars plotted
%   Both 't' and 'dat' can include columns, giving different lines
%   If the column sizes are different, excess lines are truncated
%   The value of g.minbar{n} is the minimum relative error-bar 
%   Error bars are  plotted if  relative error > 'g.minbar{n}'
%   A limit of g.cutoff{n} is used as a lower limit for log graphs.
%   An input of g.logs{n} is used to to determine whether linear or log.
%   Horizontal error-bars are allowed, with sampling + error bars combined
%   The cell array 'g.linestyle{n}' gives the line-styles that are used.
%   If there are too few line-styles, they are re-used in cyclic order
%   Licensed by Peter D. Drummond, (2022) - see License.txt 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS AND ERROR BARS

nd = size(dat);                               %%Input data size
ax = [ax(1),ax(end)];                         %%restore data axis
lg=g.legends{n};                              %%legend switch
lines=g.linestyle{n};                         %%Input linestyles
figure;                                       %%Start new graphics figure
scale = abs(g.esample{n});                    %%Set scaling factor
dat(:,:,3) = dat(:,:,3)*scale;                %%Scale sampling error
if g.esample{n} < 0                           %%If scaling factor < 0
    dat(:,:,2)=sqrt(dat(:,:,2).^2+dat(:,:,3).^2);%%Combine errors
    dat(:,:,3)=0.0*dat(:,:,3);                %%Set sampling error  zero
end                                           %%End if sampling error
nt = size(t);                                 %%Set number of time-points
et = zeros(1,nd(2));                          %%Initialize time-errors
if length(nt) == 3                            %%Check if time-errors
  if nt(3) > 2                                %%If time sampling errors
    t(:,:,3) = t(:,:,3)*abs(g.esample{n});    %%Scale sampling error 
    t(:,:,2) = t(:,:,2)+t(:,:,3);             %%Combine both error bars
  end                                         %%End if time-errors
  if nt(3) > 1                                %%Check if time step errors
    et =reshape(t(:,:,2),nt(1),nt(2));        %%Horizontal error-bar
  end                                         %%End if any step errors 
  t = reshape(t(:,:,1),nt(1),nt(2));          %%Reshape to lines+points
end                                           %%End check if time-errors 
jmax = length(lines);                         %%Number of line-styles
if nt(1) > 1 && nd(1) > 1                     %%Test for line dimension>1
    jmin = min(nt(1),nd(1));                  %%Prevent index overflow
else                                          %%Else index is set to 1
    jmin = max(nt(1),nd(1));                  %%Use maximum dimension
end                                           %%End check line dimension  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   RESHAPE INPUT DATA

for s =-1:2:1                                 %%Loop over sampling lines
 for j = 1:jmin                               %%Loop over data lines
 j1 = mod(j-1,jmax)+1;                        %%Index for line-styles
 if lines{j1} ~= 0                            %%If line-style not zero
  if j <= nd(1)                               %%If j within data range
    range = max(dat(j,:,1))-min(dat(j,:,1));
    da =reshape(dat(j,:,1),1,nd(2));          %%Line data
    max_err = max(dat(j,:,2));
    rel_err = max_err/(range+1.e-100);        %%Prevent divide by zero
    rel_et = 0;
    eu =reshape(dat(j,:,2),1,nd(2));          %%Set upper error-bar
    y =da+s*reshape(dat(j,:,3),1,nd(2));      %%Include sampling-error
    y = max(y, g.graphcutoffs{n});            %%include y cut-off
    el = max(y-eu,g.graphcutoffs{n});         %%compute lower eb location
    el = y-el;                                %%compute lower eb
  end                                         %%end if j within data range
  t1 = t;                                     %%initialise time coordinate
  sizet = size(t);
  if j <= sizet(1)
    if length(nt) == 3
      max_et =    max(et(j,:,1));
      ranget = max(t(j,:,1))-min(t(j,:,1));
      rel_et = max_et/(ranget+1.e-100);       %%Prevent divide by zero
    end
    t1 =reshape(t(j,:),1,nd(2));              %% time reshape for plotting
    et1 = reshape(et(j,:),1,nd(2));           %% time parametric errors 
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   GRAPH PLOTTING SECTION
%g.linewidth{n} =2;               %%Get n-th plot linewidth

  if  rel_et > g.minbar{n}  && nd(3) > 1          %%Parametric error bars?
    plt = errorbar(t1,y,el,eu,et1,et1,lines{j1}); %%Error-bars + sampling
  elseif  rel_err > g.minbar{n}                   %%Print error bars?
    plt = errorbar(t1,y,el,eu,lines{j1});         %%Error-bars + sampling
  else                                            %%No error bars needed!
    plt = plot(t1,y,lines{j1});                   %%2D plot
  end                                             %%End if parametric

set(plt,'linewidth',g.linewidth{n});
  hold on                                         %%Set hold for next line
 end                                              %%End if line not zero
 end                                              %%End loop over data line
 if s== -1 &&  ~isequal('',lg)                    %%If legend required 
  legend(lg,'location','best','AutoUpdate','off');%%Set legend 
 end                                              %%End if legend required 
end                                               %%End loop over sampling 
if ~g.octave                                      %%if not octave
    axis = gca;
    mint = min(min(t));
    maxt = max(max(t));
    if mint < maxt
       axis.XLim=[mint,maxt];
    end
end                                              %%end if not octave
xproperties(n,ax,g);                             %%set  graph properties
end                                              %%end xplot