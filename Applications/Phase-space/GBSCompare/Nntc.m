function C = Nntc(p)
% C = NNC(p); generates comparison count probabilities 
% for an n-fold partition.
if p.thermal ~= 1
    error('thermal inputs required for nnc comparisons');
end
partition = p.part{p.k};                         %current graph partition
maxm = 1+p.max;                                  %max boson count+1
nd = length(partition);                          %partition dimension
sz = ones(1,nd)*maxm;                            %size of probability data                                     
C = zeros([sz,1]);                               %initialize probability
p1 = zeros([1,nd]);                              %initialize numbers
M  = zeros([1,nd]);                              %initialize sizes
if iscell(partition)                             %check if cell partition
    for j =1:nd                                  %loop on partitions
        i = partition{j}(1);                     %first partition index
        p1(j)  = 1/(1+(sinh(p.sqz(i)).*p.tr).^2);   %thermal numbers
        M(j)  = length(partition{j});            %partition length
    end                                          %end loop on partitions
else                                             %if partition vector
    i = 1;                                       %first detector index
    M = partition;                               %partition lengths
    for j=1:nd                                   %loop on vector index         
        p1(j) = 1/(1+(sinh(p.sqz(i)).*p.tr).^2);    %thermal numbers
        i=i+M(j);                                %update detector index
    end                                          %end loop on vector
end                                              %end check if cell
for j =1:nd                                      %loop on partitions
    n1=maxm^(j-1);                               %array size below j
    n2=maxm^(nd-j);                              %array size above j
    C = reshape(C,n1,maxm,n2);                   %reshape log prob.
    C(:,1,:) = C(:,1,:) + M(j)*log(p1(j));       %update log prob. for j
    lp2 = log(1-p1(j));                          %get log(1-n) for j
    for m = 1:p.max                              %loop on counts
      C(:,m+1,:) = C(:,m,:)+log(1+(M(j)-1)/m)+lp2;
    end                                          %end loop on counts
end                                              %end partition loop
C = reshape(C,[sz,1]);                           %reshape log prob.
C = exp(C);                                      %probability
end                                              %end nnc
