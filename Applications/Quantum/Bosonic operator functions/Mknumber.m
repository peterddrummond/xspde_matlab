function psi = Mknumber(nvec,p)
%psi = MKNUMBER(nvec,p) makes a number state
%nvec is a vector of the numbers in each mode

nvec(end:p.modes) = nvec(end);
id = num2cell(1+nvec(1:p.modes));
psi = zeros([p.nmax,1]);
idx = sub2ind(size(psi),id{:});
psi(idx) = 1;                                  %%using sparse methods
if p.sparse                                    %%using sparse methods
  psi = reshape(psi,prod(p.nmax),1);
  if p.quantum == 2                            %%using density matrices
    psi = sparse(psi*psi');
  end
end
end