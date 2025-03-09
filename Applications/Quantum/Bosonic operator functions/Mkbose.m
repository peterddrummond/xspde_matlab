function a = Mkbose(varargin)
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
    a{j}  = kron(kron(I2,a1{j}),I1);             %%a{j} = j-th annihilation
end                                              %%end for m-mode operators
end                                              %%end mkbose

function [p,list] = mklist(varargin)
%[p,list] = MKMODES(varargin) makes p.modes, p.nmax, list

p = varargin{end};
if ~isfield(p,'nmax')                            %%reconstruct 'p.nmax'
    p.nmax = 2;
end
if ~isfield(p,'modes')                           %%reconstruct 'p.modes'
    p.modes = length(p.nmax);
end
if nargin > 1                                    %%make 'list' if not present
    list = varargin{1};
else
    list = 1:p.modes;
end
p.modes = max(p.modes,list(end));                %%increase 'p.modes' if needed
while length(p.nmax) < p.modes                   %%lengthen 'p.nmax' if needed
    p.nmax(end+1) = p.nmax(end);
end
end