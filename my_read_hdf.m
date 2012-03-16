function data = my_read_hdf(filename)

% newer QuickPIC verisons: hd5
if( findstr(filename, 'h5') )
  data = my_read_h5(filename)';
% older QuickPIC verisons: hd4
else
  fileinfo = hdfinfo(filename);
  %  Retrieve info about Scientific Data Set in example.hdf
  data_set_info = fileinfo.SDS;
  %  Check the size
  data_set_info.Dims.Size;
  % Read a subset of the data using info structure
  data = hdfread(data_set_info);
end% if

return;

