% assume particle matrix in quickpic format
%   column order: x y z x' y' p
function my_densplot(pp, hist_var1, hist_var2, do_cbar, n_sigma_axes)

if nargin<4, do_cbar = 0; end
if nargin<5, n_sigma_axes = 5; end

%
n_bins = 512;
n_colorscaling = 24 * (512/n_bins);
%n_colorscaling = 12;

pplabel(1).var = 'x [um]';
pplabel(2).var = 'y [um]';
pplabel(3).var = 'z [um]';
pplabel(4).var = 'xp [urad]';
pplabel(5).var = 'yp [urad]';
pplabel(6).var = 'p [GeV/c]';
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
hist_Zx = linspace(x_min, x_max, n_bins);
hist_Zy = linspace(E_min, E_max, n_bins);
% give 2D hist
if( size(pp,2) == 7)
% give weighted 2D hist
  [histmat, xbin, ybin] = hist2_weighted(pp(:,hist_var1), pp(:,hist_var2), hist_Zx, hist_Zy, pp(:,7) );
else
%  [histmat, xbin, ybin] = hist2(pp(:,hist_var1), pp(:,hist_var2), hist_Zx, hist_Zy);
end% if
  
h_g = pcolor(hist_Zx,hist_Zy,histmat); 
shading('flat');
dens_min = min(min(histmat));
dens_max = max(max(histmat));
caxis([dens_min dens_max/n_colorscaling]);
if(do_cbar)
  h_c = colorbar;
  set(get(h_c,'ylabel'),'String', 'density  [a.u.]');
end% if
axis([x_min x_max E_min E_max]);
%axis square tight;
xlabel(pplabel(hist_var1).var);
ylabel(pplabel(hist_var2).var);
