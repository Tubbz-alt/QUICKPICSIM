% assume particle matrix in quickpic format
%   column order: x y z x' y' p
function [hist_Z, hist_N, z_scaled] = my_histplot(pp, hist_var1, n_bins, n_sigma, do_plot)

if nargin<4, n_sigma = 5; end
if nargin<5, do_plot = 1; end

pplabel(1).var = 'x [um]';
pplabel(2).var = 'y [um]';
pplabel(3).var = 'z [um]';
pplabel(4).var = 'xp [urad]';
pplabel(5).var = 'yp [urad]';
pplabel(6).var = 'p [GeV/c]';

hz_mean = mean(pp(:,hist_var1));
hz_sigma = std(pp(:,hist_var1));
z_min = (hz_mean-n_sigma*hz_sigma);
z_max = (hz_mean+n_sigma*hz_sigma);
hist_Z = z_min:((2*n_sigma*hz_sigma+3))/n_bins:z_max;
hist_Z(1) = -inf;
hist_Z(end) = +inf; % NB: using this last bin of histc will by default have 0 counts
[hist_N bin] = histc(pp(:,hist_var1), hist_Z);
% if charge per particles column is included, we weight histogram accordingly
if( size(pp,2) == 7)
%  disp('start accum');
  hist_N = accumarray(bin, pp(:,7))'; % adjust histogram for particle weigth
  hist_N(length(hist_Z)-1) = 0; % fill up to the end if no hits
%  disp('end accum');
else
  hist_N = hist_N(1:end-1)'; % last bin if histc output will have 0 counts, ignore
end% if
% disregard outliers : cut tail points above |n_sigma*std| from dist
hist_Z = hist_Z(2:end-2); 
hist_N = hist_N(2:end-1);
% convert to mm if too large
z_scaled = 0;
if( max(abs(hist_Z)) >= 1e3)
  hist_Z = hist_Z / 1e3;
  z_min = z_min / 1e3;
  z_max = z_max / 1e3;
  pplabel(hist_var1).var(findstr(pplabel(hist_var1).var, 'u')) = 'm';  
  z_scaled = 1;
end% if
% scale y-axis to 1
hist_N = hist_N / max(hist_N);
% visualize full dist
if(do_plot)
  bar(hist_Z, hist_N);
  xlabel(pplabel(hist_var1).var);
  ylabel('count [a.u.]');
  myaxis = axis;
  axis([z_min  z_max  myaxis(3)  myaxis(4)]);
  grid on;
end% if

