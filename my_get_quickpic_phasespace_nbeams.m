function N_beams = my_get_quickpic_phasespace_nbeams(mydir)

% establish quickpic output format
qp_version_suffix = my_get_quickpic_format(mydir);
if( strfind(qp_version_suffix, '.h5') )
  beamdir = 'RAW-BEAM/';
else
  beamdir = 'PHA-BEAM/';
end% of

% writes the number of part files to temp.txt
temp_filename = [mydir beamdir 'temp.txt'];
system(['ls -lad ' mydir beamdir '*0* | wc > ' temp_filename]);
% gets number of parts
fid = fopen(temp_filename, 'r');
wc_output = str2num(fgets(fid));
N_beams = wc_output(1);
if(N_beams < 1)
  warning(['EA: no beam dirs was found.  Stopping execution.']);
  stop;
end% 
fclose(fid);

return;
