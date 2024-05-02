function a=xdst1(a,n) 
%% a=XDST1(a,n) gives a DST-I+ transform on index n
%% This is a DST-I, leaving out the two end-points
% Definition: N = input dimension, for k=2,..N-1
% a'_k = v*\sum_{n=2}^{N-1} a_n sin[pi*(n-1)(k-1)/(N-1)]
% where: v = sqrt(2/(N-1)), which makes it self-inverting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = size(a);
sx = [prod(sz(1:n-1)),sz(n),prod(sz(n+1:end))];
a  = reshape(a,sx);
n2 = 2*(sx(2)-1);
a1 = zeros(sx(1),n2,sx(3));
a1(:,2:sx(2)-1,:)  = a(:,2:sx(2)-1,:);
a1(:,sx(2)+1:n2,:) = -flip(a(:,2:sx(2)-1,:),2);
a1 = fft(a1,[],2);
a(:,2:sx(2)-1,:) = a1(:,2:sx(2)-1,:)*sqrt(2/(sx(2)-1))/(-2*1i);
a(:,1,:) = 0;
a(:,sx(2),:) = 0;
a  = reshape(a,sz);
end