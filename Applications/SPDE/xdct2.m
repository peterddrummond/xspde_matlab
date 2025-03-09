function ct=xdct2(a, tn) %% F=XSDCTN_2(a) gives a DCT-II transform on index tn
% Definition (0-based indices):
% X_k = [ \sum_{n=0}^{N-1} x_n cos[pi/N (n+1/2) k] ] / norm
% norm = sqrt(2/N)
% FFT implementation (0-based indices):
% X_k = [ real( FFT[y_n] exp(-i pi/(2N) k) ) ] / norm
% y_n = (x_0, ..., x_{N-1}, 0, ..., 0) where length(y_n) is 2*N
% The index N-1 is fixed to be zero. Hence, only the elements
% x_0,..,x_{N-2} are considered for the FFT. The number of elements N' is
% reduced to N-1.
    if ~isreal(a)           %transform real and imaginary parts seperately
    ct = xdct2(real(a),tn) + 1j*xdct2(imag(a),tn);
        return;
    end
    sz = size(a);
    sx = [prod(sz(1:tn-1)),sz(tn),prod(sz(tn+1:end))];   
    a  = reshape(a,sx);         %flatten dims before and after tn                             
    N  = size(a,2);             %length of dimension to be transformed
    M  = N-1;                   %reduced number of elements in dimension tn
    M2 = 2*M;                   %expand dimension to be transformed to M2; 
    x  = zeros(prod(sz(1:tn-1)),M2,prod(sz(tn+1:end))); %to be FFT'd    
    x(:,1:M,:) = a(:,1:N-1,:);  %copy elements       
    y = fft(x,[],2);            %carry out FFT
    postfact = exp(-1j*pi/(2*M)*(0:M-1));  %exponential post-FFT factor
    y(:,1:M,:) = y(:,1:M,:) .* postfact;   %multiply with post-FFT factor
    ct = zeros(sx);             %output tensor         
    ct(:,1:N-1,:) = real(y(:,1:M,:));  %truncate result, take real part
    ct = ct * sqrt(2/M);      %normalize
    ct = reshape(ct,sz);        %reshape to original size
end