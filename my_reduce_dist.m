% reduce # of particles by factor N_red, randomely
function pp = my_reduce_dist(pp, N_red, seed);

if nargin<3, seed = rng; end

rng(seed);

% randomize before reduce
I_rand = my_shuffle(1:size(pp,1));
pp = pp(I_rand, :);
pp = pp (1:N_red:end, :);
