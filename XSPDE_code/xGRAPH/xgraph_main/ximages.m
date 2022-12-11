function [] =  ximages(n,datan,nx,x,xlab,olab,ax,g)
%   XIMAGE_PLOT() makes 3d transverse 'movie' images at fixed axis 1
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

ax=ax(2:4);                                    %%select plotted axes
g.images{n}   = min(nx(2),g.images{n});        %%image numbers
for i = 0:g.images{n}-1                        %%Loop over images to plot
    if g.images{n} == 1                        %%only one image wanted
        np = nx(2);                            %%Print last image 
    else                                       %%several images wanted
        np = (nx(2)-1)/(g.images{n}-1);        %%image spacing in time
        np = round(1+i*np);                    %%round to next time step   
    end                                        %%End if images == 
    im = reshape(datan(1,np,:,:,1),nx(3),nx(4)); %%Image at time t(np)
    if ~isempty(g.headers{n})                        %%if full xheader wanted
        imtitle = [g.headers{n},', ',xlab{1},sprintf(' = %0.3f',x{1}(np))];
    else                                       %%end if xheaders
        imtitle= ' ';                          %%make image header blank
    end
    xpr(1,g,'Image plot in graph %d, %s = %d\n',n,xlab{1},x{1}(np));
    xplot3D(x{2},x{3},im',xlab{2},xlab{3},olab,imtitle,n,ax,g);
end                                            %%End images loop
end                                            %%End image function