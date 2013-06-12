%
% E. Adli, June 6, 2013
%    First version!
%
% MAPS: Y
function outvecX = my_2d_mapping(X, Y, invecY)
for jx=1:length(invecY),
  [n,i] = min(abs(invecY(jx)-Y));
  outvecX(jx) = X(i);
end% for


