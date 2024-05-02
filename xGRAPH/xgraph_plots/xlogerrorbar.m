function  xlogerrorbar(varargin)
%   XLOGERRORBAR(X,Y,E,min,..) plots the graph of vector X vs. vector Y with
%   a logarithmic Y-axis, using positive bars about the data points, ie:
%   plots Y plus or minus E 
%   the bars are such that the lower bar y-l is bounded below  by min, 
%   and the ratio of (y-l)/y is also bounded below by min. 
%
%   XLOGERRORBAR has a similar syntax to the original Matlab ERRORBAR.
%   The only difference is that it includes a minimum argument input.
%   History: adapted from ERRORBARLOG, Matlab Central 2021/01/05: v1.00

error(nargchk(3,inf,nargin));
e=varargin{3};
min=varargin{4};
y=max(min,varargin{2});

% displays the upper and lower error bars at y +/- e
% ymin = y-e is bounded by min

u = e;
ymin = max(min,y-e);
ymin = max(y*min,ymin);
l = y - ymin;
h = errorbar(varargin{1:2},l,u,varargin{5:end});
set(gca,'YScale','log'); % set the Y axis in log coordinates.
end