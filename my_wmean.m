function mu = my_wmean(pp, pp_var);
mu = sum(pp(:,pp_var) .* pp(:,7)) / sum(pp(:,7));
   