function a = mkbose(varargin)
%a = MKBOSE((list,pow,) p) makes annihilation operators
%One input makes all a_k operators to the first power.
%Two inputs makes a list of a_k operators to the first power.
%Three inputs makes a specified list of operators to a specified power.

pow = 1;                                         %%initial matrix power
[p,list] = mklist(varargin{:});                  %%make mode number, list
if nargin == 3                                   %%if 3 input arguments
    pow = varargin{2};                           %%reset the power
end                                              %%end if 3 input arguments
a1 = cell(p.modes); a = a1;                      %%initialize the cells
for j = list                                     %%loop on 1-mode operators
    ax = sparse(diag(sqrt(1:p.nmax(j)-1),1));    %%make 1-mode operators
    a1{j} = ax^pow;                              %%make nonlinear operators
end                                              %%end for 1-mode operators
for j = list                                     %%loop on m-mode operators
    I1  = 1; I2 = 1;                             %%Initialize
    for j1 = 1:j-1                               %%loop on lower identities
        I1  = kron(I1,speye(p.nmax(j1)));        %%I1 = j-1 mode identity
    end                                          %%end loop for identities
    for j2 = j+1:p.modes                         %%loop on upper identities
        I2  = kron(I2,speye(p.nmax(j2)));        %%I2 = m-j-1 mode identity
    end                                          %%end loop for identities
    a{j}  = kron(kron(I1,a1{j}),I2);             %%a{j} = j-th annihilation
end                                              %%end for m-mode operators
end                                              %%end mkbose