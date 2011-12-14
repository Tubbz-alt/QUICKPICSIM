function data = my_read_hdf(filename)

fileinfo = hdfinfo(filename);
%  Retrieve info about Scientific Data Set in example.hdf
data_set_info = fileinfo.SDS;
%  Check the size
data_set_info.Dims.Size;
% Read a subset of the data using info structure
data = hdfread(data_set_info);

return;

