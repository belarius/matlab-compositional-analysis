function y = clr(x)
%CLR Performs the Centered Log Ratio transform from compositional geometry into Euclidian geometry.
%   Detailed explanation goes here
%
% written by:
% Greg Jensen
% greg.guichard.jensen@gmail.com

y = log(x./repmat(geomean(x,2),1,size(x,2)));

end

%REFERENCES
%
% Aitchison, J. (1986). The statistical analysis of compositional data.
%     Chapman & Hall, Ltd.