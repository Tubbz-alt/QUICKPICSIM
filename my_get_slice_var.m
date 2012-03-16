%   slice_var1: slice along this var
%   slice_var2: calc' quantity for this var
%   slice_z_fixed: use input slice_z  (typical use: calculate for first timestep, and stick to it)
function [slice_mean,slice_sigma, slice_z, slice_N_z] = my_get_slice_var(pp, slice_var1, slice_var2, slice_z_input, N_slices, h_Nstd)

slice_z_fixed = 1;
if nargin<4, slice_z_fixed = 0; slice_z_input=0; end
if nargin<5, N_slices = 21; end
if nargin<6, h_Nstd = 3; end

% slice up beam in slice var 1
[Y, I] = sort(pp(:,slice_var1));
pp_sort = pp(I, :);
if(~slice_z_fixed)
  mu_z = mean(pp_sort(:,slice_var1));
  sigma_z = std(pp_sort(:,slice_var1));
  z_min = mu_z - h_Nstd*sigma_z;
  z_max = mu_z + h_Nstd*sigma_z;
  slice_z = [z_min:(z_max-z_min)/(N_slices):z_max];  % 2: end of 1st slice, 3: end of 2nd slice etc
else
  slice_z = slice_z_input;
end% if slice_z_fixed

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

% calc' quantities per slice 
%n_slices(2)
slice_mean(1) = mean(pp(1:n_slices(2), slice_var2));
slice_sigma(1) = std(pp(1:n_slices(2), slice_var2));
for n=2:length(n_slices),
  slice_mean(n) = mean(pp(n_slices(n-1)+1:n_slices(n), slice_var2));
  slice_sigma(n) = std(pp(n_slices(n-1)+1:n_slices(n), slice_var2));
end% for

% set NaN to 0  (will be weighted against (0) charge)
slice_mean( isnan(slice_mean) ) = 0;
slice_sigma( isnan(slice_mean) ) = 0;
