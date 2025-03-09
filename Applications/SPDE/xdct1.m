function a = xdct1(a,n) 
%% a=XDCT1(a,n) gives a DCT-I transform on index n
% Definition: N = input dimension, for k=2,..N-1
% a'_k = v*[\sum_{n=2}^{N-1} a_n cos[pi*(n-1)(k-1)/(N-1)]
%         +1/2( a_1 - (-1)^n*a_N)]
% where: v = sqrt(2/(N-1)), which makes it self-inverting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isreal(a)           %transform real and imaginary parts seperately
        a = xdct1(real(a),n) + 1j*xdct1(imag(a),n);
        return;
end
sz = size(a);
sx = [prod(sz(1:n-1)),sz(n),prod(sz(n+1:end))];
a  = reshape(a,sx);
n2 = 2*(sx(2)-1);
a1 = zeros(sx(1),n2,sx(3));
a1(:,1:sx(2),:) = a;
a1(:,sx(2)+1:n2,:) = flip(a(:,2:sx(2)-1,:),2);
a1 = fft(a1,[],2);
a  = reshape(a1(:,1:sx(2),:),sz)*sqrt(1/(2*(sx(2)-1)));
end