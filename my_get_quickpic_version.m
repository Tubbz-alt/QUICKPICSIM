% extract QuickPIC version output format from file
function data = my_get_quickpic_version(dirname)
% decide quickpic output version
temp_filename = ['/tmp/temp.txt'];
system(['ls -la ' datadir ' >| ' temp_filename]);
fid = fopen(temp_filename, 'r');
realfilename_all = fgetl(fid)
fclose(fid);
stop
