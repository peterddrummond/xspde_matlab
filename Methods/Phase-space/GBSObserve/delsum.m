function C = delsum(a,p)
% C = delv(a,p) generates the sum of u,v variances for qsim
% follows conventions from the entanglement sction

C = delv2(a,p)+delu2(a,p);
end                                              %end du function