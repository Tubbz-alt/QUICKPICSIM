function [my_std, my_mu, my_M2] = my_std(x, y);
% calc std deviation for a histogram, assuming lin interpolation
my_mu = 0;
my_M2 = 0;
N_split = 1000;
xi = min(x):(max(x)-min(x))/N_split:max(x);
yi = interp1(x,y,xi);
ea = 0;
plot(x,y,'o',xi,yi, 'x')
for(n=2:length(xi)),
  my_mu = my_mu + xi(n-1)*yi(n-1)*(xi(n)-xi(n-1));
  my_M2 = my_M2 + xi(n-1)^2*yi(n-1)*(xi(n)-xi(n-1));
end%  for
my_mu = my_mu / (max(xi)-min(xi));
my_M2 = my_M2 / (max(xi)-min(xi));
my_std = sqrt(my_M2 - my_mu^2);
