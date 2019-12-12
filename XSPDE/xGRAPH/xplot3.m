function [] =  xplot3(x1,x2,x3,lab1,lab2,lab3,imhead,imagetype)
%   XPLOT3(x1,x2,x3,lab1,lab2,lab3,imhead,imagetype)
%   Makes 3D plots depending on the imagetype
%   Here x1,x2,x3 are the data in the three dimensions
%   Note that x3 is an array
%   lab1,lab2,lab3 are labels in each of the 3 dimensions
%   imhead is a header string
%   imagetype labels the type: 1=surface, 2=greyscale, 3=contour
%   Licensed by Peter D. Drummond, (2015) - see License.txt 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
figure;                                        %%Start new graphics figure
    switch imagetype                           %%Select image type
        case 1                                 %%first image type
            surf(x1,x2,x3);                    %%plot 3d surface vs x,y
        case 2                                 %%Second image type
            contourf(x1,x2,x3);                %%plot grayplot vs x,y 
        case 3                                 %%Third image type
            contour(x1,x2,x3);                 %%plot 3d contour vs x,y
        case 4                                 %%Fourth image type
            pcolor(x1,x2,x3);                  %%plot pseudocolor vs x,y
    end                                        %%End select image type
    xheader(imhead,lab1,lab2,lab3);            %%plot title and labels
end