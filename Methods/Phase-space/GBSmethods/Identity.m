function [U] = Identity(p)
%[U] = IDENTITY(M)
%Returns an identity matrix, size M*M

if isnumeric(p)  && p <= 0
    U =  'identity';
    return;
end
U = eye(p.modes);
end