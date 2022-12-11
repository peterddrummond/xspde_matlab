function [] =  xplot3D(x1,x2,x3,lab1,lab2,lab3,imhead,n,ax,g)
%   XPLOT3D(x1,x2,x3,lab1,lab2,lab3,imhead,n,ax,g)
%   Makes 3D plots depending on the imagetype
%   Here x1,x2,x3 are the data in the three dimensions
%   Note that x3 is an array in general
%   lab1,lab2,lab3 are labels in each of the 3 dimensions
%   imhead is a header string
%   imagetypes: 1=surface, 2=greyscale, 3=contour, 4=pseudocolor
%   ax is the list of axes plotted in 3d, for calculating plot limits
%   Note that ax = in.dimension + 1 describes the observable!
%   plot properties are set using the xproperties function
%   Licensed by Peter D. Drummond, (2022) - see License.txt 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax = [ax(1:2),ax(end)];                        %%restore data axis
x3 = max(x3, g.graphcutoffs{n});               %%include z cut-off
figure;                                        %%Start new graphics figure
    switch g.imagetype{n}                      %%Select 3d image type
        case 1                                 %%First image type
            surf(x1,x2,x3);                    %%plot 3d surface vs x,y
        case 2                                 %%Second image type
            contourf(x1,x2,x3);                %%plot grayplot vs x,y 
        case 3                                 %%Third image type
            contour(x1,x2,x3);                 %%plot 3d contour vs x,y
        case 4                                 %%Fourth image type
            pcolor(x1,x2,x3);                  %%plot pseudocolor vs x,y
    end                                        %%End select image type
    xheader(imhead,lab1,lab2,lab3);            %%plot title and labels
    xproperties(n,ax,g);                       %%set plot properties
end                                            %%End function