function [] =  xtransverseplot(n,datan,nx,x,xlab,olab,g)
%   xtransverseplot() makes 2d transverse plots at fixed time
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

imhead= ' ';
g.transverse{n} = min(nx(3),g.transverse{n});  %%transverse numbers
nd = size(datan);
for i = 0:g.transverse{n}-1                    %%Loop over transverse plots
    if g.transverse{n} == 1                    %%only one plot wanted
        np = nx(3);                            %%Print last plot
    else                                       %%several plots wanted
        np = (nx(3)-1)/(g.transverse{n}-1);    %%plot spacing in time
        np = round(1+i*np);                    %%round to next time step   
    end                                        %%End if images == 1
    xp = x{2};                                 %%x-coordinate set to x[2}
    if ~isempty(g.headers{n})                  %%if full header wanted
        imhead = [olab,', ',xlab{1},sprintf(' = %0.3f',x{1}(np))];
    end                                        %%end if headers
    data = reshape(datan(:,:,np,:),[nd(1),nd(2),nd(4)]);  
    xplot(xp,data,g.esample{n},g.minbar{n},g.lines{n});
    xheader(imhead,xlab{2},olab,'');           %%title
end                                            %%End transverse plot loop
end                                            %%end Transverse function
