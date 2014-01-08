function ad = aitchdist(x,y)
%AITCHDIST Computes the distance between coordinates in compositional space.
%   Detailed explanation goes here
%
% written by:
% Greg Jensen
% Columbia University (Psychology)
% greg.guichard.jensen@gmail.com

[m,n] = size(x);
if nargin < 2
    y = ones(m,n).*(1/n);
else
    y = log(y./repmat(geomean(y,2),1,n));
end
x = log(x./repmat(geomean(x,2),1,n));
ad = (sum((x-y).^2,2)).^0.5;

end

