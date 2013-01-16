function data = my_read_h5(filename)

char_slash = findstr(filename, '/');
char_dot = findstr(filename, '.');

dataname = filename(char_slash(end):char_dot(end)-1);
dataname( dataname == '_' ) = '-'; % for some reason QuickPIC has slightly different filename than dataname

%fileinfo = h5info(filename);
% Read a subset of the data using info structure

%data = h5read(filename, dataname);
data = hdf5read(filename, dataname); % OLD LAPTOP MATLAB
return;

