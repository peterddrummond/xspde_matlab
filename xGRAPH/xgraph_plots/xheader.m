function [] =  xheader(h1,x1,x2,x3)
%   XHEADER(h1,x1,x2,x3) makes graph headers and labels axes
%   Here h1 is a header string, and x1,x2,x3 are axis label strings
%   Licensed by Peter D. Drummond, (2021) - see License.txt 

title(h1);                                     %%Set title
xlabel(x1);                                    %%Set x-axis label
ylabel(x2);                                    %%Set y-axis label
zlabel(x3);                                    %%Set z-axis label
end                                            %%end xheader