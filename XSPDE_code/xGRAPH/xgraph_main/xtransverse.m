function [] =  xtransverse(n,datan,param,nx,x,xlab,olab,ax,g)
%   xtransverse() makes 2d transverse graphs at fixed time
%   parametric plots are output if required.
%   First data dimension is the line index, last dimension is the errors
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

g.transverse{n} = min(nx(2),g.transverse{n});    %%transverse indices
nd = size(datan);
for i = 0:g.transverse{n}-1                      %%Loop over plots
    if g.transverse{n} == 1                      %%only one plot wanted
        np = nx(2);                              %%Print last plot
    else                                         %%several plots wanted
        np = (nx(2)-1)/(g.transverse{n}-1);      %%plot spacing in time
        np = round(1+i*np);                      %%round to next time step   
    end                                          %%End if images == 1
    xcoord = x{2};                               %%x-coordinate set to x[2}
    if ~isempty(g.headers{n})                    %%if full xheader wanted
        imtitle = [g.headers{n},', ',xlab{1},sprintf(' = %0.3f',x{1}(np))];
    else                                         %%end if xheaders
        imtitle= ' ';                            %%make image header blank
    end
    data = reshape(datan(:,np,:,:),[nd(1),nd(3),nd(4)]);
    xlabel = xlab{2};
    if (g.parametric{n}(2)==2) && ~isequal(param,0)%%Parametric x coord 
       xlabel = g.olabels{g.parametric{n}(1)};
       xcoord = reshape(param(:,np,:,:),[nd(1),nd(3),nd(4)]);
    end
    xpr(1,g,'2D plot in graph %d, %s = %d\n',n,xlab{1},x{1}(np));
    xplot2D(xcoord,data,n,ax,g);
    xheader(imtitle,xlabel,olab,'');             %%title
end                                              %%end plot loop
end                                              %%end xtransverse function