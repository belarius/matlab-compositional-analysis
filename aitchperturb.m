function z = aitchperturb(x,y)
%AITCHPERTURB Performs a compositional perturbation of x in terms of y. 
%   Formally, z = x?y
%   Note that x?y == y?x
%   Note that this assumes a scalar constant of 1.
%
% written by:
% Greg Jensen
% Columbia University (Psychology)
% greg.guichard.jensen@gmail.com

z = x.*y;
z = z./sum(z);

end

