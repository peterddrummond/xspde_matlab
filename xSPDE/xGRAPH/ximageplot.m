function [] =  ximageplot(n,datan,nx,x,xlab,olab,g)
%   XIMAGE_PLOT() makes 3d transverse 'movie' images at fixed axis 1
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

imhead= ' ';                                   %%make image header blank
g.images{n}   = min(nx(2),g.images{n});        %%image numbers
for i = 0:g.images{n}-1                        %%Loop over images to plot
    if g.images{n} == 1                        %%only one image wanted
        np = nx(2);                            %%Print last image 
    else                                       %%several images wanted
        np = (nx(2)-1)/(g.images{n}-1);        %%image spacing in time
        np = round(1+i*np);                    %%round to next time step   
    end                                        %%End if images == 
    im = reshape(datan(1,np,:,:,1),nx(3),nx(4)); %%Image at time t(np)
    if ~isempty(g.headers{n})                  %%if full xheader wanted
        imhead = [olab,', ',xlab{1},sprintf(' = %0.3f',x{1}(np))];
    end                                        %%end if xheaders
    xplot3(x{2},x{3},im',xlab{2},xlab{3},olab,imhead,g.imagetype{n});
end                                            %%End images loop
end                                            %%End image function