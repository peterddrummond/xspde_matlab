function [] =  xplot(t,datan,esample,minbar,lines)
%   XPLOT(t,datan,ebar,minbar)
%   Makes 2D plots with error bars and sampling errors against axis 't'
%   Here 'datan' is the n-th graph data, with sampling error and error bars
%   Note: 't' is a dummy label for axis-1, it could be anything
%   Input data has three dimensions: lines, points and errors
%   Input t has up to three dimensions, but at least one
%   The value of 'esample' indicates how sampling error is graphed
%   If esample = 0, no sampling error lines are plotted, just the mean
%   If esample = -n, n*sigma sampling errors are included in errorbars
%   If esample = n, upper + lower n*sigma sampling error lines are plotted
%   Both 't' and 'datan' can include columns, giving different lines
%   If the column sizes are different, excess lines are truncated
%   The value of minbar is the minimum relative error-bar 
%   Error bars are  plotted if  relative error > 'minbar'
%   Horizontal error-bars are allowed, with sampling + error bars combined
%   The cell array 'lines' gives the line-styles that are used.
%   If there are too few line-styles, they are re-used in cyclic order
%   Licensed by Peter D. Drummond, (2015) - see License.txt 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PARAMETERS AND ERROR BARS

nx = size(datan);                             %%Input data size
figure;                                       %%Start new graphics figure
datan(:,:,3) = datan(:,:,3)*abs(esample);     %%Scale sampling error
if esample < 0                                %%If sampling error <0
      datan(:,:,2)=datan(:,:,2)+datan(:,:,3); %%Combine both error bars
      datan(:,:,3)=0.0*datan(:,:,3);          %%Set sampling error  zero
end                                           %%End if sampling error
nt = size(t);                                 %%Set number of time-points
et = zeros(1,nx(2));                          %%Initialize time-errors
if length(nt) == 3                            %%Check if time-errors
  if nt(3) > 2                                %%If sampling errors
    t(:,:,3) = t(:,:,3)*abs(esample);         %%Scale sampling error 
    t(:,:,2) = t(:,:,2)+t(:,:,3);             %%Combine both error bars
  end                                         %%End if time-errors
  if nt(3) > 1                                %%Check if any step errors
    et =reshape(t(:,:,2),nt(1),nt(2));        %%Horizontal error-bar
  end                                         %%End if any step errors 
  t = reshape(t(:,:,1),nt(1),nt(2));          %%Reshape to lines+points
end                                           %%End check if time-errors 
jmax = length(lines);                         %%Number of line-styles
if nt(1) > 1 && nx(1) > 1                     %%Test for line dimension>1
    jmin = min(nt(1),nx(1));                  %%Prevent index overflow
else                                          %%Else index is set to 1
    jmin = max(nt(1),nx(1));                  %%Use maximum dimension
end                                           %%End check line dimension                                   %%End test for line dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   RESHAPE INPUT DATA???

for j = 1:jmin                                %%Loop over line index
 j1 = mod(j-1,jmax)+1;                        %%Index for line-styles
 if lines{j1} ~= 0                            %%If line-style not zero
  if j <= nx(1)
    max_error = max(datan(j,:,2));
    max_range = max(datan(j,:,1))-min(datan(j,:,1));
    rel_error = max_error/(max_range+1.e-100);%%Prevent divide by zero
    rel_et = 0;
    da =reshape(datan(j,:,1),1,nx(2));        %%Transverse data 
    eb =reshape(datan(j,:,2),1,nx(2));        %%Transverse error-bar
    se =reshape(datan(j,:,3),1,nx(2));        %%Sampling-error
  end
  t1 = t;
  sizet = size(t);
  if j <= sizet(1)
    if length(nt) == 3
      max_et =    max(et(j,:,1));
      max_ranget = max(t(j,:,1))-min(t(j,:,1));
      rel_et = max_et/(max_ranget+1.e-100);   %%Prevent divide by zero
    end
    t1 =reshape(t(j,:),1,nx(2));              %% time
    et1 = reshape(et(j,:),1,nx(2));           %% time parametric errors 
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   GRAPH PLOTTING SECTION

  if  rel_et > minbar                         %%Parametric error bars?
    errorbar(t1,da+se,eb,eb,et1,et1,lines{j1}); %%Error-bars + sampling
    hold on;                                  %%Allows combined figures
    errorbar(t1,da-se,eb,eb,et1,et1,lines{j1}); %%Error-bars - sampling
  elseif  rel_error > minbar                  %%Print error bars?
    errorbar(t1,da+se,eb,lines{j1});          %%Error-bars + sampling
    hold on;                                  %%Allows combined figures
    errorbar(t1,da-se,eb,lines{j1});          %%Error-bars - sampling
  else                                        %%No error bars!
    plot(t1,da+se,lines{j1});                 %%2D plot, black, upper sd
    hold on;                                  %%Allows combined figures
    plot(t1,da-se,lines{j1});                 %%2D plot, black, lower sd
  end                                         %%End if errorchecks
 end
end
mint = min(min(t));
maxt = max(max(t));
if mint>maxt
     xlim([mint,maxt]);                       %%Set plot limits
end
end                                           %%end xplot