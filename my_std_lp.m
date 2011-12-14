function [sig, mu] = my_std_lp(x, y);
% calc std deviation for a histogramdist

N_split = 10000;
xi = min(x):(max(x)-min(x))/N_split:max(x);
yi = interp1(x,y,xi);
%xi = x;
%yi = y;
%plot(xi,yi, '-x');
mu = sum(xi.*yi) / sum(yi)
sig = sqrt( 1 / (sum(yi)-1) *  sum( yi.*(xi - mu).^2) )
