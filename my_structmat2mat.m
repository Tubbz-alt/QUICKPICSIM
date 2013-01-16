% Usage :
%     function AxB = my_calc_cross(AX, AY, AZ, BX, BY, BZ);
% 1) A x B of a matrix
% 2) |A x B| 
% 3) (A x B) / |B|^2
%
%% Changelog :
%% E. Adli, May 3, 2012 
%%   First version!
function mymat = my_structmat2mat(mystructmat, myindex);

if( nargin < 2 )
  myindex = 1;
end% if

  tempvar = struct2cell(mystructmat);
  tempvar2 = reshape(tempvar(myindex,:,:), size(tempvar(myindex,:,:), 2), size(tempvar(myindex,:,:),3));
  mymat = cell2mat(tempvar2);
