% assume particle matrix in quickpic format
%   column order: x y z x' y' p
function my_densplot(pp, hist_var1, hist_var2, do_cbar, n_sigma_axes)

if nargin<4, do_cbar = 0; end
if nargin<5, n_sigma_axes = 5; end

n_bins = 2^ceil(log(sqrt(size(pp, 1))) / log(2) );
n_bins = 2^ceil(log(sqrt(size(pp, 1))) / log(2) / 1.15); % better dude
n_colorscaling = 24 * (512/n_bins);

pplabel(1).var = 'x [um]';
pplabel(2).var = 'y [um]';
pplabel(3).var = 'z [um]';
pplabel(4).var = 'xp [urad]';
pplabel(5).var = 'yp [urad]';
pplabel(6).var = 'p [GeV/c]';
pplabel(8).var = 'p_{meas} [GeV/c]';
%pplabel(9).var = '\Delta p [GeV/c]';
pplabel(9).var = '(p_{meas} - p_{real}) / p_{real} [%]';
x_min = mean(pp(:,hist_var1))-n_sigma_axes*std(pp(:,hist_var1));
x_max = mean(pp(:,hist_var1))+n_sigma_axes*std(pp(:,hist_var1));
E_min = mean(pp(:,hist_var2))-n_sigma_axes*std(pp(:,hist_var2));
E_max = mean(pp(:,hist_var2))+n_sigma_axes*std(pp(:,hist_var2));
if( max(abs(x_min), abs(x_max)) >= 1e3)
  x_min = x_min / 1e3;
  x_max = x_max / 1e3;
  pp(:,hist_var1) = pp(:,hist_var1) / 1e3;
  pplabel(hist_var1).var(findstr(pplabel(hist_var1).var, 'u')) = 'm';  
end% if  
if( max(abs(E_min), abs(E_max)) >= 1e3)
  E_min = E_min / 1e3;
  E_max = E_max / 1e3;
  pp(:,hist_var2) = pp(:,hist_var2) / 1e3;
  pplabel(hist_var2).var(findstr(pplabel(hist_var2).var, 'u')) = 'm';  
end% if  
hist_Zx = linspace(x_min*0.99, x_max*1.01, n_bins);
hist_Zy = linspace(E_min*0.99, E_max*1.01, n_bins);
% give 2D hist
if( size(pp,2) == 7)
% give weighted 2D hist
  [histmat, xbin, ybin] = hist2_weighted(pp(:,hist_var1), pp(:,hist_var2), hist_Zx, hist_Zy, pp(:,7) );
else
  [histmat, xbin, ybin] = hist2(pp(:,hist_var1), pp(:,hist_var2), hist_Zx, hist_Zy);
end% if

% convert to log
%histmat = log10(histmat);

% custom color map
addpath('/Users/eadli/Dropbox/SLAC/E200/E200_shared/E200_scripts/tools/');
cmap = custom_cmap();
colormap(cmap.wbgyr);

h_g = pcolor(hist_Zx,hist_Zy,histmat); 
shading('flat');
%dens_min = min(min(histmat))
dens_min = 0;
dens_max_all = sort(max(histmat));
dens_max = dens_max_all(round(end*95/100));
%caxis([dens_min dens_max/n_colorscaling]); % better without ....
caxis([dens_min dens_max]); % better without ....
caxis([dens_min dens_max]); % better without ....
if(do_cbar)
  h_c = colorbar;
  set(get(h_c,'ylabel'),'String', 'density  [a.u.]');
end% if
axis([x_min x_max E_min E_max]);
%axis square tight;
xlabel(pplabel(hist_var1).var);
ylabel(pplabel(hist_var2).var);

%colormap('jet');
