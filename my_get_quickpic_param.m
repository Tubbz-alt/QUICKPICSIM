function myval = my_get_quickpic_param(filename, myparam, n_occ_req)

if( nargin<3 )
  n_occ_req = 1; % default, pick first occurence
end% if

n_occ = 0;

fid = fopen(filename, 'r');
row = 0;
n_row = 0;
match_found = 0;
while( row ~= -1 )
  n_row = n_row + 1;
  row = fgets(fid);
  if( row ~= -1)
    n_col = findstr(row, myparam);
    if (n_col>0) 
      % parse to "="
      while( (n_col < length(row)) & (row(n_col) ~= '=') )
        n_col = n_col + 1;
      end% while
      n_col = n_col + 1;
      % parse to end if value
      myval = [];
      while( (n_col < length(row)) & (row(n_col) ~= ',') & (row(n_col) ~= '/')  )
        myval = [myval row(n_col)];
        n_col = n_col + 1;
      end% while
      if(exist('myval') & ischar(myval) & ~isempty(str2num(myval)) )
          match_found = 1; % pick first real match (later will be ignored)
          n_occ = n_occ + 1;
      end% if
    end% if
  end%if
  if( match_found && (n_occ == n_occ_req) )
    break;
  end% if
end% while
fclose(fid);


if(exist('myval'))
  if(ischar(myval) &&  (n_occ == n_occ_req))
    myval = str2num(myval);
  else
    myval = -1;
  end% if
else
  myval = -1;
end% if

return;
