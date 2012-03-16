% extract QuickPIC version output format from file
function qp_version_suffix = my_get_quickpic_format(dirname)

% decide quickpic output version
temp_filename = ['/tmp/temp.txt'];

system(['ls -la ' dirname ' >| ' temp_filename]);

% gets number of parts
fid = fopen(temp_filename, 'r');
ls_output = 0;
qp_version_suffix = [];
while( ls_output ~= -1 )
ls_output = fgets(fid);
if( ls_output ~= -1)
  if( findstr('RAW-BEAM', ls_output) )
    qp_version_suffix = '.h5';
  end% if
  if( findstr('PHA-BEAM', ls_output) )
      qp_version_suffix = '.hdf';
  end% if
end% if
end% while
fclose(fid);

