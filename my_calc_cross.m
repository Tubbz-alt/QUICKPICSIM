% Usage :
%     function AxB = my_calc_cross(AX, AY, AZ, BX, BY, BZ);
% 1) A x B of a matrix
% 2) |A x B| 
% 3) (A x B) / |B|^2
%
%% Changelog :
%% E. Adli, May 3, 2012 
%%   First version!
function AxB = my_calc_cross(AX, AY, AZ, BX, BY, BZ);

for ix=1:size(AX,1);
  for jx=1:size(AY,2);
    A = [AX(ix, jx) AY(ix, jx) AZ(ix, jx)];
    B = [BX(ix, jx) BY(ix, jx) BZ(ix, jx)];
    AxB(ix, jx).vec = cross(A, B);
    AxB(ix, jx).X = AxB(ix, jx).vec(1);
    AxB(ix, jx).Y = AxB(ix, jx).vec(2);
    AxB(ix, jx).Z = AxB(ix, jx).vec(3);
    AxB(ix, jx).mag = norm(AxB(ix, jx).vec, 2);
    AxB(ix, jx).vec_v = AxB(ix, jx).vec / norm(B, 2)^2;
    AxB(ix, jx).mag_v = norm(AxB(ix, jx).vec_v, 2);
  end% for
end% for


% NOTE: using matlabs vector cross was not faster
if(0)
for ix=1:size(AX,1);
  for jx=1:size(AY,2);
    A(jx, :) = [AX(ix, jx) AY(ix, jx) AZ(ix, jx)];
    B(jx, :) = [BX(ix, jx) BY(ix, jx) BZ(ix, jx)];
    AxB(jx, :).vec = cross(A, B, 2);
  end% for
end% for
%    AxB(jx).mag = norm(AxB(ix, jx).vec, 2);
%    AxB(jx).vec_v = AxB(ix, jx).vec / norm(B, 2)^2;
%    AxB(jx).mag_v = norm(AxB(ix, jx).vec_v, 2);
%  end% for
%end% for
end% if
