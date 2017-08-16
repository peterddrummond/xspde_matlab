function [] =  xplot(t,datan,esample,minbar,lines)
%   XPLOT(t,datan,ebar,minbar)
%   Makes 2D plots with error bars and sampling errors against time t
%   Here 'datan' is the n-th graph data, with sampling error and error bars
%   The value of 'esample' indicates how sampling error is graphed
%   If esample = 0, no sampling error lines are plotted, just the mean
%   If esample = -n, nsigma sampling errors are included in the errorbars
%   If esample = n, upper and lower nsigma sampling error lines are plotted
%   The value of minbar is the minimum relative error-bar 
%   Error bars are  plotted if  relative error > 'minbar'
%   The cell array 'lines' gives the line-styles that are used.
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

nx = size(datan);
figure;                                       %%Start new graphics figure
datan(:,3,:) = datan(:,3,:)*abs(esample);     %%Scale sampling error 
if esample < 0
      datan(:,2,:)=datan(:,2,:)+datan(:,3,:);  %%Combine both error bars
      datan(:,3,:)=0.0*datan(:,3,:);           %%Set sampling error to zero
end                                            %%End if sampling error            

for j = 1:nx(1)
  max_error = max(datan(j,2,:));
  max_range = max(datan(j,1,:))-min(datan(j,1,:));
  rel_error = max_error/(max_range+1.e-100);   %%Prevent divide by zero   
  da =reshape(datan(j,1,:),1,nx(3));           %%Transverse data 
  eb =reshape(datan(j,2,:),1,nx(3));           %%Transverse error-bar
  se =reshape(datan(j,3,:),1,nx(3));           %% sampling-error                       
  if  rel_error > minbar                       %%Print error bars?
    errorbar(t,da+se,eb,lines{j});             %%Error-bars + sampling
    hold on;                                   %%Allows combined figures
    errorbar(t,da-se,eb,lines{j});             %%Error-bars - sampling
  else                                         %%No error bars!
    plot(t,da+se,lines{j});                    %%2D plot, black, upper sd
    hold on;                                   %%Allows combined figures
    plot(t,da-se,lines{j});                    %%2D plot, black, lower sd
  end                                          %%End if errorchecks
end
xlim([t(1),t(length(t))]);                     %%Set plot limits
end                                            %%end xplot2