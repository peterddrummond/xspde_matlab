function [] =  xplot2(t,da,eb,se,ebar,minbar)
%   XPLOT(t,da,eb,se,checks,minbar)
%   Makes 2D plots with error bars and sampling errors
%   Here 'da' is data, 'se' is sampling error, 'eb' is an error bar
%   Next, 'checks' indicates if error bars should be printed
%   Limits lim1, lim2 are the upper and lower x-axis limits
%   Error bars are only printed if relative size exceeds 'minbar'
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

figure;                                        %%Start new graphics figure
relative_error = max(eb)/(max(da)-min(da));
if  ebar  && (relative_error > minbar)         %%Are there error bars?
    errorbar(t,da+se,eb,'k');                  %%Error-bars + sampling
    hold on;                                   %%Allows combined figures
    errorbar(t,da-se,eb,'k');                  %%Error-bars - sampling
else                                           %%No error bars!
    plot(t,da+se,'k');                         %%2D plot, black, upper sd
    hold on;                                   %%Allows combined figures
    plot(t,da-se,'k');                         %%2D plot, black, lower sd
end                                            %%End if errorchecks
xlim([t(1),t(length(t))]);                     %%Set plot limits
end                                            %%end start xplot