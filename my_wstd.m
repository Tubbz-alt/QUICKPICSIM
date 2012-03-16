function sig = my_wstd(pp, pp_var);
mu = my_wmean(pp, pp_var);
sig = sqrt( sum((pp(:,pp_var) - mu).^2.*pp(:,7) )  / (sum(pp(:,7))) );
   