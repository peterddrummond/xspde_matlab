function ct=xdct3(a, tn) 
%% F=XDCT3(a, tn) gives a DCT-III transform on index tn
% Definition (0-based indices):
% X_k = [ 1/2 x_0 + \sum_{n=1}^{N-1} x_n cos[pi/N (k+1/2) n] ] / norm
% norm = sqrt(2/N)
% FFT implementation (0-based indices):
% X_k  = [ real( FFT[y_n] ) + 1/2 x_0 ] / norm
% y_n  = (0, xx_1,..., xx_{N-1}, 0, ..., 0) where length(y_n) is 2*N
% xx_n = x_n exp(-i pi/(2 N) n)
% The index N-1 is fixed to be zero. Hence, only the elements
% x_0,..,x_{N-2} are considered for the FFT. The number of elements N' is
% reduced to N-1.
    if ~isreal(a)           %transform real and imaginary parts seperately
        ct = xdct3(real(a),tn) + 1j*xdct3(imag(a),tn);
        return;
    end
    sz = size(a);
    sx = [prod(sz(1:tn-1)),sz(tn),prod(sz(tn+1:end))];   
    a  = reshape(a,sx);      %flatten dims before and after tn                             
    N  = size(a,2);          %length of dimension to be transformed
    M  = N-1;                %reduced number of elements in dimension tn
    M2 = 2*M;                %expand dimension to be transformed to M2;     
    x  = zeros(prod(sz(1:tn-1)),M2,prod(sz(tn+1:end))); %to be FFT'd    
    x(:,1:M,:) = a(:,1:N-1,:);  %copy elements    
    x(:,1,:) = 0;            %remove 1st element from FFT
    prefact = reshape(exp(-1j*pi/(2*M)* (0:M-1)),[1,M,1]);  %exponential pre-FFT factor
    x(:,1:M,:) = x(:,1:M,:) .* prefact;                     %multiply with pre-FFT factor
    y = fft(x,[],2);         %carry out FFT
    ct = zeros(sx);          %output tensor
    ct(:,1:N-1,:) = real(y(:,1:M,:)); %truncate result, take real part
    ct(:,1:N-1,:) = ct(:,1:N-1,:) + 1/2*a(:,1,:);   %add 1st element
    ct = ct * sqrt(2/M);     %normalize
    ct = reshape(ct,sz);     %reshape to original size
end
