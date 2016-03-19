function [] =  xtransverseplot(n,datan,nx,x,xlab,olab,g)
%   xtransverseplot() makes 2d transverse plots at fixed time
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

imhead= ' ';
g.transverse{n} = min(nx(1),g.transverse{n});  %%transverse numbers
for i = 0:g.transverse{n}-1                    %%Loop over transverse plots
    if g.transverse{n} == 1                    %%only one plot wanted
        np = nx(1);                            %%Print last plot
    else                                       %%several plots wanted
        np = (nx(1)-1)/(g.transverse{n}-1);    %%plot spacing in time
        np = round(1+i*np);                    %%round to next time step   
    end                                        %%End if images == 1
    xp = x{2};                                 %%x-coordinate set to x[2}
    if ~isempty(g.headers{n})                  %%if full header wanted
        imhead = [olab,', ',xlab{1},sprintf(' = %0.3f',x{1}(np))];
    end                                        %%end if headers
    da_n =reshape(datan(1,np,:),1,nx(2));      %%Transverse data 
    eb_n =reshape(datan(2,np,:),1,nx(2));      %%Transverse error-bar    
    se_n =reshape(datan(3,np,:),1,nx(2));      %%Transverse sampling-error     
    x2plot(xp,da_n,eb_n,se_n,g.ebar{n},g.minbar{n});
    xheader(imhead,xlab{2},olab,'');            %%title
end                                            %%End transverse plot loop
end                                            %%end Transverse function
