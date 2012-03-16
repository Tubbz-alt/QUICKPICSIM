%   slice_var1: slice along this var
%   slice_var2: calc' quantity for this var
%   slice_z_fixed: use input slice_z  (typical use: calculate for first timestep, and stick to it)
function [emnx,emny,betax,alphax,betay,alphay] = my_get_slice_twiss(pp, slice_z_input)

slice_z_fixed = 1;
N_slices = 21; 
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

% calc' twiss quantities per slice 
n_range_1 = 1:n_slices(2);
pp_placet = [pp(n_range_1,6)/1e1 pp(n_range_1,1) pp(n_range_1,2) pp(n_range_1,3) pp(n_range_1,4) pp(n_range_1,5)];
TW = my_calc_twiss(pp_placet);
eps_N = my_calc_emittance(pp_placet); % not yet macro particle
emnx(1) = eps_N(1);
emny(1) = eps_N(2);
betax(1) = TW(1,1);
alphax(1) = -TW(1,2);
betay(1) = TW(3,3);
alphay(1) = -TW(3,4); 

for n=2:length(n_slices),
n_range = (n_slices(n-1)+1):n_slices(n);
pp_placet = [pp(n_range,6)/1e1 pp(n_range,1) pp(n_range,2) pp(n_range,3) pp(n_range,4) pp(n_range,5)];
TW = my_calc_twiss(pp_placet);
eps_N = my_calc_emittance(pp_placet); % not yet macro particle
emnx(n) = eps_N(1);
emny(n) = eps_N(2);
betax(n) = TW(1,1);
alphax(n) = -TW(1,2);
betay(n) = TW(3,3);
alphay(n) = -TW(3,4); 
end% for

% set NaN to 0  (will be weighted against (0) charge)
emnx( isnan(emnx) ) = 0;
emny( isnan(emny) ) = 0;
betax( isnan(betax) ) = 0;
alphax( isnan(alphax) ) = 0;
betay( isnan(betay) ) = 0;
alphay( isnan(alphay) ) = 0;

