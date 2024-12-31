function [CF] = xqchooseftnc(P, p)
%  CF = xxqchooseftnc(P,p) gives weighted multinomial coefficients. 
%  This treats the case where P is an Md-component partition cell vector.
%  Each component, P{n}, is a vector of indices, allowing arbitrary indices
%  All indices listed in each P{n} are then combined into one partition
%  which has a total size of M(n) = length(P{n}) for the n-th partition.
%  If the first dimension of p is a scalar it assumes p is a constant
%  and returns the product of successive coefficients of x^n(i) in the 
%  binomial expansion of (xp+(1-p))^M(n) for each of the Md dimensions. 
%  If first dimension of p is an Mn-vector it assumes p is inhomogeneous, 
%  and returns successive coefficients of x^n(i) in prod_j(xp_j+(1-p_j))
%  This selects the number of times that p=1 in each partition
%  If p has a second dimension, it averages the result over this.
%  Results are accurate to 15 decimals at most, due to IEEE round-off

%INITIALIZE VECTORS AND CONSTANTS

Md = length(P);                             % Get dimension of partition
Mn = 0;                                     %? Initialize  max index
P1 = zeros(1, Md);                          % Create vector of sizes
for L = 1:Md                                % Loop over the partitions
	Mn = max(Mn,max(P{L}));                 % Get  max index in P{L}
	P1(L) = length(P{L}) + 1;               % Generate vector of sizes
end                                         % End loop over the partitions
sz = size(p);                               % Get input size of p matrix
if sz(1) == 1                               % Check if first size is 1
    p(1:Mn,:) = p(1,:);                     % If 1, expand p to Mn terms 
end                                         % End check if first size is 1
theta = 2*pi./P1;                           % vector of Fourier angles
Mlp = [P1, sz(2)];                          % Extended P1 vector of sizes
PI = ones(Mlp);                             % Initialize data array

%INDEX OVER THE INPUT DATA
        
for L = 1:Md                                % Loop over data partitions 
    Mx(1) = prod(Mlp(1:L-1));               % Lower M dimensions of output
    Mx(2) = Mlp(L);                         % Middle M dimensions of output
    Mx(3) = prod(Mlp(L+1:Md));              % Upper M dimensions of output
    Mx(4) = sz(2);                          % Ensemble dimension of output
    PI = reshape(PI, Mx);                   % Reshape data array to MX
    pr = p(P{L},:);                         % Restrict input to a partition 
    for k = 1:Mlp(L)                        % Loop on k inside partition
        expk = exp(-1i*(k-1)*theta(L));     % Exponential factor for k     
        cpk = (1-pr) + pr.*expk;            % Get projector * exponential
        probs(1,1,1,:) = prod(cpk,1);       % Get k-product over partition
        PI(:,k,:,:) = PI(:,k,:,:).*probs;   % Replace jth col. with prod.
    end                                     % End loop inside partition                    
end                                         % End loop over partitions 
CF = reshape(PI, Mlp);                      % Reshape array to full size
%CF = mean(PI, Md+1);                        % Perform averaging

%INVERSE FOURIER TRANSFORM OVER EACH PARTITION

for n1 = 1:Md                              % Loop over partitions 
    CF = ifft(CF, P1(n1), n1);             % Inverse FFT  of partition n1
end                                        % End loop over partitions 
CF = real(CF);                             % Return only real data
end                                        % End function