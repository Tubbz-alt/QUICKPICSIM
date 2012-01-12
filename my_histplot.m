% assume particle matrix in quickpic format
%   column order: x y z x' y' p
function [XXz_total, NNz_total] = my_histplot(pp, hist_var1, n_bins, n_sigma, do_plot)

if nargin<4, n_sigma = 5; end
if nargin<5, do_plot = 1; end

pplabel(1).var = 'x [um]';
pplabel(2).var = 'xp [urad]';
pplabel(3).var = 'y [um]';
pplabel(4).var = 'yp [urad]';
pplabel(5).var = 'z [um]';
pplabel(6).var = 'p [GeV/c]';

hz_mean = mean(pp(:,hist_var1));
hz_sigma = std(pp(:,hist_var1));
z_min = (hz_mean-n_sigma*hz_sigma);
z_max = (hz_mean+n_sigma*hz_sigma);
hist_z = z_min:(2*n_sigma*hz_sigma)/n_bins:z_max;
[NNz_total, XXz_total] = hist(pp(:,hist_var1), hist_z);
% disregard outliers
NNz_total(1) = 0;
NNz_total(end) = 0;
% visualize full dist
if(do_plot)
  bar(XXz_total, NNz_total);
  xlabel(pplabel(hist_var1).var);
  ylabel('count [a.u.]');
  myaxis = axis;
  axis([z_min  z_max  myaxis(3)  myaxis(4)]);
  grid on;
end% if

