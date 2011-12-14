function qp_BEAM = my_get_quickpic_phasespace(mydir, n_beam_str, n_timestep_str)

qp_BEAM = [];

% writes the first line to get starting number (it seems to not always be 0 for beam 2)
temp_filename = [mydir 'PHA-BEAM/' n_beam_str '/' 'temp.txt'];
system(['ls -la ' mydir 'PHA-BEAM/' n_beam_str '/PHA-BEAM-' '*-' n_beam_str '_' n_timestep_str '* > ' temp_filename]);
% gets number of parts
fid = fopen(temp_filename, 'r');
ls_output = 0;
while( ls_output ~= -1 )
ls_output = fgets(fid);
  if( ls_output ~= -1)
    n_str = findstr('PHA-BEAM-', ls_output);
    if(n_str == [])
      warning(['EA: part file not found for ' 'PHA-BEAM-****' '-' n_beam_str '_' n_timestep_str '.hdf' '.  Stopping execution.']);
      stop;
    end% 
    part_filename = ls_output(n_str:n_str+24);
    myfile = [mydir 'PHA-BEAM/' n_beam_str '/'  part_filename];
    qp_BEAM_part = double(my_read_hdf(myfile));
    qp_BEAM = [qp_BEAM; qp_BEAM_part];
  end% if
end% while
fclose(fid);

return;
