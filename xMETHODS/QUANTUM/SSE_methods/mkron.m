function psi = mkron(M,varargin)
%   psi = mkron(M,varargin) Computes multiple kronecker products
%   The input is M1 vectors v1,v2,..
%   The output is psi(i_1,i_2,..i_M) = v1(i_1)*v2(i_2)...
%   If there are M1 < M input vectors, the last one is repeated
%   If there are M1 > M input vectors, the last one(s) are ignored
%   xSDE functions licensed by Peter D. Drummond, (2023) - see License

M1 = min(nargin-1,M);
fields = zeros(1,M);
tsize  = cell(1,M);
for i = 1:M
  if i <= M1
    fields(i) = length(varargin{i});
  else
    fields(i) = fields(i-1);
    varargin{i}=varargin{i-1};
  end
  tsize{i} = [ones(1,i-1),fields(i),ones(1,M-i+1)];
end 
psi = ones([fields,1]);
for i = 1:M
   psi = psi.*reshape(varargin{i},tsize{i});
end
end  