% assume particle matrix in quickpic format
%   column order: x y z x' y' p
function my_densplot(pp, hist_var1, hist_var2, do_cbar, n_sigma_axes)

if nargin<4, do_cbar = 0; end
if nargin<5, n_sigma_axes = 5; end

pplabel(1).var = 'x [um]';
pplabel(2).var = 'y [um]';
pplabel(3).var = 'z [um]';
pplabel(4).var = 'xp [urad]';
pplabel(5).var = 'yp [urad]';
pplabel(6).var = 'p [GeV/c]';
x_min = mean(pp(:,hist_var1))-n_sigma_axes*std(pp(:,hist_var1));
x_max = mean(pp(:,hist_var1))+n_sigma_axes*std(pp(:,hist_var1));
E_min = mean(pp(:,hist_var2))-3*std(pp(:,hist_var2));
E_max = mean(pp(:,hist_var2))+3*std(pp(:,hist_var2));
xedges = linspace(x_min, x_max, 512);
yedges = linspace(E_min, E_max, 512);

histmat = hist2(pp(:,hist_var1), pp(:,hist_var2), xedges, yedges);
h_g = pcolor(xedges,yedges,histmat); 
shading('flat');
dens_min = min(min(histmat));
dens_max = max(max(histmat));
caxis([dens_min dens_max/12]);
if(do_cbar)
  h_c = colorbar;
  set(get(h_c,'ylabel'),'String', 'density  [a.u.]');
end% if
axis([x_min x_max E_min E_max]);
%axis square tight;
xlabel(pplabel(hist_var1).var);
ylabel(pplabel(hist_var2).var);
