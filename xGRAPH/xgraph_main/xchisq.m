function   xchisqplot(n,d,g) 
%   XCHISQPLOT(n,d,g) plots chi-squared data  
%   for a single xGRAPH plot function n, with probability data
%   Requires comparison data in the inputs, ie last index = 4
%   Input: graph number 'n', data 'd', graph parameters 'g'.
%   Output: graphs of chi-squared data 
%   Sums over all the data points in a sampled distribution
%   Can cmake one-sided or two-sided chi-square comparisons
%   Results are graphed against any space-time dimensions
%   Called by: xgraph
%   Needs:     xmultigraph
%   Licensed by Peter D. Drummond, (2021) - see License.txt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SUM OVER CHI-SQUARE FOR PLOTS
% 
dimension = 0;                                   %%initialize the dimension
if isfield(g,'dimension')                        %%check if dimension field
    dimension = g.dimensions ;                   %%store the dimension
end                                              %%end check if field
if isequal(dimension,0)                          %%if a plot not possible?
    return                                       %%return to calling code
end                                              %%end if plot
sz         = size(d);                            %% data size
axsize     = length(sz) - 2;                     %%total axes available
nsumd      = (axsize-dimension+2):(axsize+1);    %%indices for summing
nplotd     = 1:(axsize-dimension+1);             %%indices for plotting
sz1        = prod(sz(nplotd));                   %%total plot size
sz2        = prod(sz(nsumd));                    %%total summing size
d          = reshape(d,[sz1,sz2,6]);             %%reshape count data
if g.scale{n} > 0                                %%standard chi-squared
    scale  = g.scale{n};                         %%count scaling
    var    = d(:,:,4)*scale;                     %%variance from counts
    cutoff = g.mincount;                         %%count cutoff
else                                             %%two-sided chi-squared
    scale  = 1;                                  %%no count scaling
    var    = d(:,:,3).^2 +d(:,:,6).^2;           %%variance from data
    cutoff = g.cutoffs{n};                       %%probability cutoff
end                                              %%end if scale
count      = d(:,:,1)*scale;                     %%scale simulated counts
exp        = d(:,:,4)*scale;                     %%scale expected counts
nc = (count>cutoff)&(exp>cutoff);                %%nonzero counts
valid      = sum(nc(:,:,1),2);                   %%get number of points
if g.chisqplot{n} > 1                            %%If plot 2 required?
    valid  = sqrt(2*valid-1);                    %%modified comparison 
end                                              %%End if plot 2 required
valid(:,:,2) = 0;                                %%set error bars to zero
valid(:,:,3) = 0;                                %%sampling errors to zero
valid = reshape(valid,[sz(nplotd),3]);           %%reshape counts data
chisqu = nc.*((count-exp).^2)./(var+1.e-100);    %%get chi-squared data
chisqu = sum(chisqu,2);                          %%sum over probabilities
chisqu(:,:,2) = 0;                               %%set error bars to zero
chisqu(:,:,3) = 0;                               %%sampling errors to zero
chisqu = reshape(chisqu,[sz(nplotd),3]);         %%reshape chi-square data
if g.chisqplot{n} == 1                           %%First plot type?
    g.olabels{n}= '\chi^2, k';                   %%make label
    chisqu = [chisqu;valid];                     %%include count data
else                                             %%Second polt type?
    g.olabels{n}='(2\chi^2)^{1/2}-(2k-1)^{1/2}'; %%make label
    chisqu = sqrt(2*chisqu)-valid;               %%subtract count data
end                                              %%End if plot type
xmultigraph(n,chisqu,0,g);                       %%make chi-square plots
end                                              %%end function call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    END XCHISQPLOT FUNCTION