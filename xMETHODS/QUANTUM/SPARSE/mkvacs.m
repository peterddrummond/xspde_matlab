function [vac] = mkvacs(p)
%[vac] = MKVACS(p) makes a vacuum state
if p.sparse
    d = prod(p.nmax);
    if p.quantum == 1                            %%if using wave-functions
        vac = zeros(d,1);                        %%zero vector
    elseif p.quantum == 2                        %%if using density matrices
        vac = sparse(d,d);                       %%zero sparse matrix
    end                                          %%end if wave-functions
    vac(1,1) = 1;                                %%occupy the ground state
else
    vac = zeros(p.nmax,1);                       %%zero vector  
    vac(ones(1+length(p.nmax),1)) = 1;           %%occupy the ground state
end
end                                              %%end mkvacs