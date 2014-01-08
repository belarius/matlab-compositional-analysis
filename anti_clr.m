function x = anti_clr(c_x)
%ANTI_CLR Reverses clr(x), resulting in a set of compositions.
%   Detailed explanation goes here
%
% written by:
% Greg Jensen
% greg.guichard.jensen@gmail.com

[~,n] = size(c_x);
x = exp(c_x);
x = x./repmat(sum(x,2),1,n);

end

