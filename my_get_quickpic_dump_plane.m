% extract QuickPIC dump plane from file
function qp_dump_plane = my_get_quickpic_dump_plane(dirname)

% decide quickpic output version
temp_filename = ['/tmp/temp.txt'];

system(['ls -la ' dirname ' > ' temp_filename]);

% gets number of parts
fid = fopen(temp_filename, 'r');
ls_output = 0;
qp_dump_plane = [];
while( ls_output ~= -1 )
ls_output = fgets(fid);
if( ls_output ~= -1)
  if( findstr('QEB-XZ', ls_output) )
    qp_dump_plane = 'XZ';
  end% if
  if( findstr('QEB-YZ', ls_output) )
      qp_dump_plane = 'YZ';
  end% if
end% if
end% while
fclose(fid);

 