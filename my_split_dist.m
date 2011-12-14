% split bunch dist into two dists at z = z_split
function [pp_1, pp_2, N_pp_1, N_pp_2] = my_split_dist(pp, z_split)

[Y, I] = sort(pp(:,5));
pp_sort = pp(I, :);

n = 1;
while(pp_sort(n,5) < z_split)
  n=n+1;
end% while

n_split = n;

pp_1 = pp_sort(1:n_split, :);
pp_2 = pp_sort(n_split+1:end, :);

N_pp_1 = length(1:n_split);
N_pp_2 = length(n_split+1:length(pp));


