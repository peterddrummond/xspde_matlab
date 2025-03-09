function [CF] = xqchooseftn(P, p)
%  CF = xqchooseftn(P,p) gives weighted multinomial coefficients. 
%  This treats the case where P is an Md-component partition vector
%  Inputs p are divided into successive partitions of length P(i).
%  If the first dimension of p is a scalar it assumes p is a constant
%  and returns the product of successive coefficients of x^n(i) in the 
%  multinomial expansion of (xp+(1-p))^P(i) for each of the Md dimensions. 
%  If first dimension of p is an Mn-vector it assumes p is inhomogeneous, 
%  and returns successive coefficients of x^n(i) in prod_j(xp_j+(1-p_j))
%  This selects the number of times that p=1 in each partition
%  If p has a second dimension, it averages the result over this.
%  Results are accurate to 15 decimals at most, due to IEEE round-off

%INITIALIZE VECTORS AND CONSTANTS

Md = length(P);                             % Get dimension of partition
Mn = sum(P);                                % Get number of points in P
sz = size(p);                               % Get input size 
if sz(1) == 1                               % Check if first size is 1
    p(1:Mn,:) = p(1,:);                     % Expand to Mn terms 
end                                         % End check if first size is 1
P1 = P(1,:) + 1;                            % Generates (P+1)-vector
theta = 2*pi./P1;                           % get angles in Fourier space 
Mlp = [P1, sz(2)];                          % Extends P1 vector
PI = ones(Mlp);                             % Initialize data array

%INDEX OVER THE INPUT DATA

ind = 1;                                    % Initialize start index           
for L = 1:Md                                % Loop over data partitions 
    Mc = P(L);                              % Current partition dimension
    Mx(1) = prod(Mlp(1:L-1));               % Lower M dimensions of output
    Mx(2) = Mlp(L);                         % Middle M dimensions of output
    Mx(3) = prod(Mlp(L+1:Md));              % Upper M dimensions of output
    Mx(4) = sz(2);                          % Ensemble dimension of output
    PI = reshape(PI, Mx);                   % Reshape data array to MX
    for j = 1:Mlp(L)                        % Loop inside each partition
        expk = exp(-1i*(j-1)*theta(L));     % Exponential term    
        pr = p(ind:ind+Mc-1,:);             % Restrict input to a partition 
        cpk = (1-pr) + pr.*expk;            % Get projector * exponential
        probs(1,1,1,:) = prod(cpk,1);       % Get product over factors
        PI(:,j,:,:) = PI(:,j,:,:).*probs;   % Replace jth col. with prod.
    end                                     % End loop inside partition 
    ind = ind + Mc;                         % Next start index                      
end                                         % End loop over partitions 
PI = reshape(PI, Mlp);                      % Reshape array to full size
%CF = mean(PI, Md+1);                        % Perform averaging
CF = PI;

%INVERSE FOURIER TRANSFORM OVER EACH PARTITION

for n1 = 1:Md                              % Loop over partitions 
    CF = ifft(CF, P1(n1), n1);             % Inverse FFT  of partition n1
end                                        % End loop over partitions 
CF = real(CF);                             % Return only real data
end                                        % End function