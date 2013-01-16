% slope > 0: find neg to pos zero crossing
% slope < 0: find pos to neg zero crossing
function n = my_find_zero_cross(f, slope)
n_zc = -1;
if(slope > 0)
  myinv = 1;
else
  myinv = -1;
end% if

for n=1:(length(f)-1);
  if( (myinv*sign(f(n)) < 0 ) && ( myinv*sign(f(n+1)) > 0) )
    n_zc = n;
    break;
  end% if
end% for

   
