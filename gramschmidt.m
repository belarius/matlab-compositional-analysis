function [Q,R] = gramschmidt(A)
%GRAMSCHMIDT Performs the Gram-Schmidt operation to produce an orthonormal basis.
%   Detailed explanation goes here
%
% written by:
% Greg Jensen
% greg.guichard.jensen@gmail.com

[~,n] = size(A);
% compute QR using Gram-Schmidt
R = zeros(n-1,n);
Q = zeros(n+1,n);
for j = 1:n
   v = A(:,j);
   for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
   end
   R(j,j) = norm(v);
   Q(:,j) = v/R(j,j);
end

end

