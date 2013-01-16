%   slice_var1: slice along this var
%   slice_var2: calc' quantity for this var
%   slice_z_fixed: use input slice_z  (typical use: calculate for first timestep, and stick to it)
function [slice_gx, slice_gxp, slice_gy, slice_gyp] = my_get_slice_gauss(pp, slice_z_input, N_slices)

if nargin<3, N_slices = 51; end
slice_z_fixed = 1;
h_Nstd = 3; 

% rely on previous slicing (assumed same beam, sorted in z)
slice_var1 = 3; % assume slice in z for twiss
[Y, I] = sort(pp(:,slice_var1));
pp_sort = pp(I, :);
slice_z = slice_z_input;
% slice up beam in n
n_slices(1) = 1;
slice_N_z = zeros(1,N_slices);
n_pp_count = 1;
for n=1:N_slices,
  if( n_pp_count <= length(pp_sort) )
    while( (pp_sort(n_pp_count,slice_var1) < slice_z(n+1)) && ( n_pp_count < length(pp_sort) ) )
      n_slices(n) = n_pp_count;
      slice_N_z(n) = slice_N_z(n) + 1;
      n_pp_count = n_pp_count + 1;
    end% while
  else
    n=N_slices + 1; % end scan, no more particles
  end% if
end% for

slice_N_z = slice_N_z(1:length(n_slices));

slice_gx = zeros(length(n_slices), 1);
slice_gxp = zeros(length(n_slices), 1);
slice_gy = zeros(length(n_slices), 1);
slice_gyp = zeros(length(n_slices), 1);

% calc' gauss fit per slice 
n_range_1 = 1:n_slices(2);
if( length(n_range_1) > 0 )
  pp_slice = pp_sort(n_range_1, :);
  [sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp_slice, 0);
  slice_gx(1) = gauss_sigx;
  slice_gy(1) = gauss_sigy;
  slice_gxp(1) = gauss_sigxp;
  slice_gyp(1) = gauss_sigyp;
  end
  


for n=2:length(n_slices),
  n_range = (n_slices(n-1)+1):n_slices(n);
if( length(n_range) > 0 )
  pp_slice = pp_sort(n_range, :);
  [sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp_slice, 0);
  slice_gx(n) = gauss_sigx;
  slice_gy(n) = gauss_sigy;
  slice_gxp(n) = gauss_sigxp;
  slice_gyp(n) = gauss_sigyp;
end% if
end% for

% set NaN to 0  (will be weighted against (0) charge)
slice_gx( isnan(slice_gx) ) = 0;
slice_gy( isnan(slice_gy) ) = 0;
slice_gxp( isnan(slice_gxp) ) = 0;
slice_gyp( isnan(slice_gyp) ) = 0;

% extrapolate with zero to length of input slice.z
if(slice_z_fixed)
  if( length(n_slices) < length(slice_z_input) )
    slice_gx( (length(n_slices)+1):length(slice_z_input)-1) = 0;
    slice_gy( (length(n_slices)+1):length(slice_z_input)-1) = 0;
    slice_gxp( (length(n_slices)+1):length(slice_z_input)-1) = 0;
    slice_gyp( (length(n_slices)+1):length(slice_z_input)-1) = 0;
  end% if
end% if
