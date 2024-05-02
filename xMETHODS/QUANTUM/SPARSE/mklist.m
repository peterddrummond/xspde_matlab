function [p,list] = mklist(varargin)
%[p,list] = MKMODES(varargin) makes p.modes, p.nmax, list

p = varargin{end};
if ~isfield(p,'nmax')                            %%reconstruct 'p.nmax'
    p.nmax = 2;
end
if ~isfield(p,'modes')                           %%reconstruct 'p.modes'
    p.modes = length(p.nmax);
end
if nargin > 1                                    %%make 'list' if not present
    list = varargin{1};
else
    list = 1:p.modes;
end
p.modes = max(p.modes,list(end));                %%increase 'p.modes' if needed
while length(p.nmax) < p.modes                   %%lengthen 'p.nmax' if needed
    p.nmax(end+1) = p.nmax(end);
end
end