function c_x = closure(x)
%CLOSURE Converts the coordinates in X into relative proportions.
%
% written by:
% Greg Jensen
% greg.guichard.jensen@gmail.com

c_x = x./repmat(sum(x,2),1,size(x,2));

end

%REFERENCES
%
% Aitchison, J. (1986). The statistical analysis of compositional data.
%     Chapman & Hall, Ltd.
% Jensen, G. (Submitted). The compositional analysis of choice: Behavior in
%     the simplex.