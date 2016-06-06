function [] =  ximageplot(n,datan,nx,x,xlab,olab,g)
%   XIMAGE_PLOT() makes 3d transverse images at fixed axis 1
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

imhead= ' ';                                   %%make image header blank
g.images{n}   = min(nx(3),g.images{n});        %%image numbers
for i = 0:g.images{n}-1                        %%Loop over images to plot
    if g.images{n} == 1                        %%only one image wanted
        np = nx(3);                            %%Print last image 
    else                                       %%several images wanted
        np = (nx(3)-1)/(g.images{n}-1);        %%image spacing in time
        np = round(1+i*np);                    %%round to next time step   
    end                                        %%End if images == 
    im = reshape(datan(1,1,np,:,:),nx(4),nx(5)); %%Image at time t(np)
    if ~isempty(g.headers{n})                  %%if full xheader wanted
        imhead = [olab,', ',xlab{1},sprintf(' = %0.3f',x{1}(np))];
    end                                        %%end if xheaders
    figure;                                    %%Start new graphics figure
    switch g.imagetype{n}                      %%Select image type
        case 1                                 %first image type
            surf(x{2},x{3},im');               %%plot 3d surface vs x,y
        case 2                                 %%Second image type
            contourf(x{2},x{3},im');           %%plot grayplot vs x,y 
        case 3                                 %%Third image type
            contour(x{2},x{3},im',10);         %%plot 3d contour vs x,y
    end                                        %%End select image type
    xheader(imhead,xlab{2},xlab{3},olab);      %%title
end                                            %%End images loop
end                                            %%End image function