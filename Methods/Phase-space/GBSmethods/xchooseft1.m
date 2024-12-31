function [CF] = xqchooseft1(M,pr)
%  CF = xqchooseft1(M,pr) calculates weighted binomial coefficients.
%  If pr is a scalar it returns the successive coefficients of x^n in the 
%  binomial expansion of (x*pr+(1-pr))^M. If pr is an M-vector, 
%  it returns successive coefficients of x^n in prod_j(x*pr_j+(1-pr_j))
%  If pr has a second dimension, it averages the result over this.
%  Here M is the length of the vector that is used. 
%  Results are accurate to at most 15 decimals, due to round-off

%CF=zeros(1,M+1);                                 %initialize correlation
theta=2*pi/(M+1);                                %initialize phase angle
sz = size(pr);                                   %initialize input size
if sz(1) == 1                                    %if first size equals 1
        pr(1:M,:) = pr(1,:);                     %expand p to M terms
end                                              %end if first size
for k=1:M+1                                      %loop over Fourier index
    cpk=(1-pr)+pr*exp(-1i*theta*(k-1));          %Get M factors p(j,:)
    probs = prod(cpk,1);                         %Product over M factors
    %CF(k) = mean(probs,2);                       %Average over ensemble
    CF(1,k,:) = probs;
end                                              %end loop over Fourier
CF=real(ifft(CF));                               %inverse Fourier transform
end                                              %end function

%  test case: suppose pr = 0.5, then nchooseft(M,pr) ~ 2^(-M)*nchoosek(M,k)
%  where nchoosek(M,k) is the Matlab binomial coefficient function
%   CF = qchooseft(M,0.5)
%   for k=0:M
%       p2 = nchoosek(M,k)*2^(-M);
%       fprintf( 'k %d, ft %15.12e, binomial %15.12e, error %15.12e\n',...
%       k,CF(k+1),p2, abs(CF(k+1)-p2));
%   end