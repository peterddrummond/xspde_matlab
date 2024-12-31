function C = delprod(a,p)
% C = delprod(a,p) generates the product of u,v standard deviations
% follows conventions from the entanglement sction

C = delv(a,p).*delu(a,p);
end                                              %end du function