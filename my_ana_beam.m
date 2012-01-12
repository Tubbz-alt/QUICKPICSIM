%
% quickpicsim_anabeam(beam, plot_vector, n_sigma) : general beam analysis function
%
% E. Adli, Jan 11, 2012
%
% input: 6x6 beam in quickpic (EA) format
%   for elegant input: convert with for example my_ele2qp()
% output: sigmas, betas etc.,  see source
function [sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam(pp, plot_vector, n_sigma);

if nargin<2, plot_vector = [0 0 0 0]; end
if nargin<3, n_sigma = 3; end

% for histograms
n_bins = 101;

SI_c = 299792458;
% to Elegant (org) format - where everything is in [m] etc
pp = my_qp2ele(pp);
% to Elegant (EA) format
pp(:,5) = pp(:,5)*SI_c; % [m]
% ot placet format
pp_placet = [pp(:,6)*0.511/1e3 pp(:,1) pp(:,3) pp(:,5) pp(:,2) pp(:,4)];
TW = my_calc_twiss(pp_placet);
eps_N = my_calc_emittance(pp_placet);
emnx = eps_N(1);
emny = eps_N(2);
betax = TW(1,1);
alphax = -TW(1,2);
betay = TW(3,3);
alphay = -TW(3,4); 
std_E = std(pp(:,6));
muE = mean(pp(:,6));
sigz = std(pp(:,5));
mu_z = mean(pp(:,5));
pp(:,5) = pp(:,5) - mu_z;

N_plot_cols = sum(plot_vector);
n_plot = 0;

%
% z E
%
CC = corrcoef( pp(:,5)*1e6, pp(:,1)*1e6);
corzx = CC(1,2);
CC = corrcoef( pp(:,5)*1e6, pp(:,3)*1e6);
corzy = CC(1,2);
CC = corrcoef( pp(:,5)*1e6, pp(:,6)*1e6);
corzp = CC(1,2);
muE = muE*0.511/1e3;
sigEE = std_E/muE; 
if(plot_vector(1) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols,n_plot);
  my_densplot(pp * diag([1 1 1 1 1e6 0.511/1e3]), 5, 6, 0, n_sigma);
  title(['\sigma_z [um]=' num2str(sigz*1e6, '%.1f') ', \mu_E [GeV] =' num2str(muE, '%.2f') ]);
  xlabel(['z [um]' ', c_{zp}=' num2str(corzp, '%.1E') ', c_{zx}=' num2str(corzx, '%.1E')  ]);
  ylabel(['E [GeV]' ', \sigma_E / E [-] =' num2str(sigEE, '%.2E') ]);
  subplot(2,N_plot_cols, n_plot + N_plot_cols);
end% if
% hist
[fit_X, fit_Y] = my_histplot(pp * diag([1 1 1 1 1e6 1]), 5, n_bins, n_sigma, plot_vector(1));
[gauss_sigz, gauss_mu_z, gauss_A_z] = my_gaussplot(fit_X, fit_Y, plot_vector(1));
if(plot_vector(1) == 1)
  title(['\sigma_{z,gauss} [um]=' num2str(gauss_sigz*1e6, '%.1f')]);
end% if

%
% x
%
sigx = std(pp(:,1));
if(plot_vector(2) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols,n_plot);
  plot(pp(:,1)*1e3, pp(:,2)*1e3, 'x')
  my_densplot(pp * diag([1e6 1e6 1 1 1 1]), 1, 2, 0, n_sigma);
  xlabel(['x [um]' ', \epsilon_{Nx} [um] =' num2str(emnx*1e6, '%.2f') ]);
  ylabel('xp [urad]');
  title(['\sigma_x [um]=' num2str(sigx*1e6, '%.1f') ', \beta_{x}=' num2str(betax, '%.3f' ) ', \alpha_{x}=' num2str(alphax, '%.3f' )  ]);
  subplot(2,N_plot_cols, n_plot + N_plot_cols);
end% if
[fit_X, fit_Y] = my_histplot(pp * diag([1e6 1 1 1 1 1]), 1, n_bins, n_sigma, plot_vector(2));
[gauss_sigx, gauss_mu_x, gauss_A_x] = my_gaussplot(fit_X, fit_Y, plot_vector(2));
if(plot_vector(2) == 1)
  title(['\sigma_{x,gauss} [um]=' num2str(gauss_sigx*1e6, '%.1f')]);
end% if


%
% y
%
sigy = std(pp(:,3));
if(plot_vector(3) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols,n_plot);
  my_densplot(pp * diag([1 1 1e6 1e6 1 1]), 3, 4, 0, n_sigma);
  xlabel(['y [mm]' ', \epsilon_{Ny} [um] =' num2str(emny*1e6, '%.2f') ]);
  ylabel('yp [mrad]');
  title(['\sigma_y [um]=' num2str(sigy*1e6, '%.1f') ', \beta_{y}=' num2str(betay, '%.3f' ) ', \alpha_{y}=' num2str(alphay, '%.3f' )  ]);
  subplot(2,N_plot_cols, n_plot + N_plot_cols);
end% if
[fit_X, fit_Y] = my_histplot(pp * diag([1 1 1e6 1 1 1]), 3, n_bins, n_sigma, plot_vector(3));
[gauss_sigy, gauss_mu_y, gauss_A_y] = my_gaussplot(fit_X, fit_Y, plot_vector(3));
if(plot_vector(3) == 1)
  title(['\sigma_{y,gauss} [um]=' num2str(gauss_sigy*1e6, '%.1f')]);
end% if


%
% xp
%
sigxp = std(pp(:,2));
if(plot_vector(4) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols, n_plot);
end% if
[fit_X, fit_Y] = my_histplot(pp * diag([1 1e6 1 1 1 1]), 2, n_bins, n_sigma, plot_vector(4));
[gauss_sigxp, gauss_mu_xp, gauss_A_xp] = my_gaussplot(fit_X, fit_Y, plot_vector(4));
if(plot_vector(4) == 1)
  title(['\sigma_{xp,gauss} [urad]=' num2str(gauss_sigxp*1e6, '%.1f')]);
end% if


%
% yp
%
sigyp = std(pp(:,4));
if(plot_vector(4) == 1)
  subplot(2,N_plot_cols, n_plot+N_plot_cols);
end% if
[fit_X, fit_Y] = my_histplot(pp * diag([1 1 1 1e6 1 1]), 4, n_bins, n_sigma, plot_vector(4));
[gauss_sigyp, gauss_mu_yp, gauss_A_yp] = my_gaussplot(fit_X, fit_Y, plot_vector(4));
if(plot_vector(4) == 1)
  title(['\sigma_{yp,gauss} [urad]=' num2str(gauss_sigyp*1e6, '%.1f')]);
end% if






