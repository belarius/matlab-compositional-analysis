function z = aitchpower(x,alpha)
%AITCHPOWER Performs a power transformation of x in terms of alpha
%   Formally, z = a?x
%   Note that a?x ~= x?a and furthermore that x?a is an invalid statement
%   Note that this assumes a scalar constant of 1.
%
% written by:
% Greg Jensen
% greg.guichard.jensen@gmail.com

z = x.^alpha;
z = z./sum(z);

end

