function aip = aitchip(x,y)
%AITCHIP Computes the inner product of coordinates in compositional space.
%   Detailed explanation goes here
%
% written by:
% Greg Jensen
% Columbia University (Psychology)
% greg.guichard.jensen@gmail.com

[~,n] = size(x);
x = log(x./repmat(geomean(x,2),1,n));
y = log(y./repmat(geomean(y,2),1,n));
aip = sum(x.*y,2);

end

