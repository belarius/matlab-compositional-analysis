function y = ilr(x,Q,type)
%ILR Performs the Isometric Log Ratio transform from compositional geometry into Euclidian geometry.
%   *** For a vector x of length len, corresponding to coordinates in
%       Aitchison geometry, ilr(x) returns a 1×(len-1) vector of
%       coordinates in Euclidian geometry.
%   *** For a matrix x of size m×n, performs ilr(x(i,:)) for each row i,
%       returning a m×(n-1) matrix.
%   B    = Either a binary partition matrix (default) or an orthonormal basis.
%   type = A string denoting 'partition' if Q is a partition matrix B
%       or 'basis' if Q is an orthonormal basis U.
%   If that Q cannot be used to initialize a coherent orthonormal basis,
%       the Gram-Schmidt procedure is used to generate a default basis.
%
% written by:
% Greg Jensen
% greg.guichard.jensen@gmail.com

[~,n] = size(x);
if nargin < 2
    V = eye(n);
    V = V(1:n,1:n-1) - [zeros(1,n-1);V(2:n,2:n)];
    U = local_gramschmidt(V);
else
    %==SET BASIS==
    if nargin==2
        type = 'partition';
    end
    if strcmp(type,'partition')
        pos = Q==1;
        neg = Q==-1;
        r = repmat(sum(pos,2),1,n);
        s = repmat(sum(neg,2),1,n);
        U = (pos.*(sqrt(s./(r.*(r+s)))) - neg.*(sqrt(r./(s.*(r+s)))))';
    elseif strcmp(type,'basis')
        U = Q;
    end
    %==TEST BASIS==
    if round(norm(U).*1000000) ~= 1000000                                   % round() is invoked to bypass imprecision in Matlab's 
        warning('ilr:IllDefinedBasis','ilr():Orthonomal basis is not well-defined; using default basis')
        V = eye(n);
        V = V(1:n,1:n-1) - [zeros(1,n-1);V(2:n,2:n)];
        U = local_gramschmidt(V);
    end
end
y = local_clr(x)*U;

end

function y = local_clr(x)
%CLR Performs the Centered Log Ratio transform from Aitchison (i.e.
%   compositional i.e. barycentric) geometry into Euclidian geometry.

[~, n] = size(x);
y = log(x./repmat(geomean(x,2),1,n));

end

function [Q,R] = local_gramschmidt(A)
%GRAMSCHMIDT Generates a default orthonormal basis

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

%REFERENCES
%
% Aitchison, J. (1986). The statistical analysis of compositional data.
%     Chapman & Hall, Ltd.
% Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal,
%     C. (2003). Isometric logratio transformations for compositional data
%     analysis. Mathematical Geology, 35(3), 279-300.
% Egozcue, J. J., & Pawlowsky-Glahn, V. (2006). Simplicial geometry for
%     compositional data. Geological Society, London, Special Publications,
%     264(1), 145-159.
% Jensen, G. (Submitted). Compositions and their application to the
%     analysis of choice.