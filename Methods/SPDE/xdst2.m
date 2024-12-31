function ct=xdst2(a, tn) 
%% F=XDST2(a,tn) gives a DST-II transform on index tn
% Definition (0-based indices):
% X_k = [ \sum_{n=0}^{N-1} x_n sin[pi/N (n+1/2) (k+1)] ] / norm
% norm = sqrt( 2 / N ) )
% FFT implementation (0-based indices):
% X_k  = [ -imag( exp(-i pi (k+1)/(2 N)  ) FFT[y_n] ) ] / norm
% y_n  = (xx_0, ..., xx_{N-1}, 0, ..., 0) where length(y_n) is 2*N
% xx_n = x_n exp(-i pi n/N)) 
% The index 0 is fixed to be zero. Hence, only the elements
% x_1,..,x_{N-1} are considered for the FFT. The number of elements N' is
% reduced to N-1.
    if ~isreal(a)           %transform real and imaginary parts seperately
        ct = xdst2(real(a),tn) + 1j*xdst2(imag(a),tn);
        return;
    end
    sz = size(a);
    sx = [prod(sz(1:tn-1)),sz(tn),prod(sz(tn+1:end))];   
    a  = reshape(a,sx);      %flatten dims before and after tn            
    N  = size(a,2);          %length of dimension to be transformed
    M  = N-1;                %reduced number of elements in dimension tn
    M2 = 2*M;                %expand dimension to be transformed to M2;       
    x  = zeros(prod(sz(1:tn-1)),M2,prod(sz(tn+1:end))); %to be FFT'd
    x(:,1:M,:) = a(:,2:N,:); %copy elements       
    prefact  = reshape(exp(-1j*pi/M*(0:M-1))   , [1,M,1]);  %exponential pre-FFT factor
    postfact = reshape(exp(-1j*pi/(2*M)* (1:M)), [1,M,1]);  %exponential post-FFT factor
    x(:,1:M,:) = x(:,1:M,:) .* prefact;                     %multiply with pre-FFT factor
    y = fft(x,[],2);         %carry out FFT
    y(:,1:M,:) = y(:,1:M,:).*postfact;  %multiply with post-FFT factor
    ct = zeros(sx);          %output tensor
    ct(:,2:N,:) = -imag(y(:,1:M,:));   %truncate result, take imag part
    ct = ct * sqrt(2/M);     %normalize
    ct = reshape(ct,sz);     %reshape to original size
end