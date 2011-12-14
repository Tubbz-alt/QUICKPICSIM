function N_beams = my_get_quickpic_phasespace_nbeams(mydir)
% writes the number of part files to temp.txt
temp_filename = [mydir 'PHA-BEAM/temp.txt'];
system(['ls -lad ' mydir 'PHA-BEAM/*0* | wc > ' temp_filename]);
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
