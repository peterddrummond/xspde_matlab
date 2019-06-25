function [] =  xtransverseplot(n,datan,param,nx,x,xlab,olab,g)
%   xtransverseplot() makes 2d transverse plots at fixed time
%   parametric plots are output if required.
%   First dimension is the field index, last dimension is the error-bars
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

imhead= ' ';
g.transverse{n} = min(nx(2),g.transverse{n});    %%transverse plot numbers
nd = size(datan);
for i = 0:g.transverse{n}-1                      %%Loop over plots
    if g.transverse{n} == 1                      %%only one plot wanted
        np = nx(2);                              %%Print last plot
    else                                         %%several plots wanted
        np = (nx(2)-1)/(g.transverse{n}-1);      %%plot spacing in time
        np = round(1+i*np);                      %%round to next time step   
    end                                          %%End if images == 1
    xcoord = x{2};                               %%x-coordinate set to x[2}
    if ~isempty(g.headers{n})                    %%if full header wanted
        imhead = [olab,', ',xlab{1},sprintf(' = %0.3f',x{1}(np))];
    end                                          %%end if headers
    data = reshape(datan(:,np,:,:),[nd(1),nd(3),nd(4)]);
    xlabel = xlab{2};
    if g.parametric{n}(2) == 2                   %%Parametric plot, x coord 
       xlabel = g.olabels{g.parametric{n}(1)};
       xcoord = reshape(param(:,np,:,:),[nd(1),nd(3),nd(4)]);
    end
    xplot(xcoord,data,g.esample{n},g.minbar{n},g.lines{n});
    xheader(imhead,xlabel,olab,'');              %%title
end                                              %%end plot loop
end                                              %%end xtransverse function
