%
% quickpicsim_anabeam(beam, plot_vector, n_sigma) : general beam analysis function
%
% E. Adli, Jan 11,  2012
%
% input: 6x6 beam in quickpic (EA) format
%   for elegant input: convert with for example my_ele2qp()
% output: sigmas, betas etc.,  see source
function [sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp, plot_vector, n_sigma);

if nargin<2, plot_vector = [1 1 1 0 0 1]; end
if nargin<3, n_sigma = 3; end

% we expand by a weight vector of ones, if none is present (overhead in calc is small)
if(size(pp,2) == 6)
  pp = [pp ones(size(pp, 1), 1)];
end% if

if(length(plot_vector) < 9)
  plot_vector(9) = 0;
end% if

if(sum(plot_vector) > 0)
  clf;
  set(0,'defaultaxesfontsize',13);
  set(0,'defaultAxesFontName', 'Arial');
end% if

% for histograms
n_bins = 201;

SI_c = 299792458;
% temp to "placet rev" format to use some pre-made functions
pp_placet = [pp(:,6)/1e1 pp(:,1) pp(:,2) pp(:,3) pp(:,4) pp(:,5) pp(:,7)];
TW = my_calc_twiss(pp_placet);
eps_N = my_calc_emittance(pp_placet); % not yet macro particle
emnx = eps_N(1);
emny = eps_N(2);
betax = TW(1,1);
alphax = -TW(1,2);
betay = TW(3,3);
alphay = -TW(3,4); 
muE = my_wmean(pp, 6);
std_E = my_wstd(pp, 6);
sigEE = std_E/muE; 
sigz = my_wstd(pp, 3);
mu_z = my_wmean(pp, 3);
pp(:,3) = pp(:,3) - mu_z;

N_plot_cols = sum(plot_vector);
n_plot = 0;


%
% z E
%
CC = corrcoef( pp(:,3).*pp(:,7), pp(:,1).*pp(:,7));
corzx = CC(1,2);
CC = corrcoef( pp(:,3).*pp(:,7), pp(:,2).*pp(:,7));
corzy = CC(1,2);
CC = corrcoef( pp(:,3).*pp(:,7), pp(:,6).*pp(:,7));
corzp = CC(1,2);
if(plot_vector(1) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols,n_plot);
  my_densplot(pp, 3, 6, 0, n_sigma);
  subplot(2,N_plot_cols, n_plot + N_plot_cols);
end% if
% hist
[fit_X, fit_Y, x_scaled] = my_histplot(pp, 3, n_bins, n_sigma, plot_vector(1));
[gauss_sigz, gauss_mu_z, gauss_A_z] = my_gaussplot(fit_X, fit_Y, x_scaled, plot_vector(1));

%
% x
%
sigx = my_wstd(pp, 1);
if(plot_vector(2) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols,n_plot);
  my_densplot(pp, 1, 4, 0, n_sigma);
  subplot(2,N_plot_cols, n_plot + N_plot_cols);
end% if
[fit_X, fit_Y, x_scaled] = my_histplot(pp, 1, n_bins, n_sigma, plot_vector(2));
[gauss_sigx, gauss_mu_x, gauss_A_x] = my_gaussplot(fit_X, fit_Y, x_scaled, plot_vector(2));
%mean(pp(:,1))
%mean(fit_X.*fit_Y)
%(fit_X(3)-fit_X(2))

%
% y
%
sigy = my_wstd(pp, 2);
if(plot_vector(3) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols,n_plot);
  my_densplot(pp, 2, 5, 0, n_sigma);
  subplot(2,N_plot_cols, n_plot + N_plot_cols);
end% if
[fit_X, fit_Y, x_scaled] = my_histplot(pp, 2, n_bins, n_sigma, plot_vector(3));
[gauss_sigy, gauss_mu_y, gauss_A_y] = my_gaussplot(fit_X, fit_Y, x_scaled, plot_vector(3));


%
% xp
%
sigxp = my_wstd(pp, 4);
if(plot_vector(4) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols, n_plot);
end% if
[fit_X, fit_Y] = my_histplot(pp, 4, n_bins, n_sigma, plot_vector(4));
[gauss_sigxp, gauss_mu_xp, gauss_A_xp] = my_gaussplot(fit_X, fit_Y, x_scaled, plot_vector(4));


%
% yp
%
sigyp = my_wstd(pp, 5);
if(plot_vector(4) == 1)
  subplot(2,N_plot_cols, n_plot+N_plot_cols);
end% if
[fit_X, fit_Y] = my_histplot(pp, 5, n_bins, n_sigma, plot_vector(4));
[gauss_sigyp, gauss_mu_yp, gauss_A_yp] = my_gaussplot(fit_X, fit_Y, x_scaled, plot_vector(4));


%
% x,E and x,y plots
%
if(plot_vector(5) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols, n_plot);
  my_densplot(pp, 1, 6, 0, n_sigma);
%  my_densplot(pp, 4, 5, 0, n_sigma); % EA VERY TEMP
  subplot(2,N_plot_cols, n_plot+N_plot_cols);
  my_densplot(pp, 1, 2, 0, n_sigma);
end% if


%
% extra z,x and z,y plots
%
if(plot_vector(6) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols, n_plot);
  my_densplot(pp, 3, 1, 0, n_sigma);
  subplot(2,N_plot_cols, n_plot+N_plot_cols);
  my_densplot(pp, 3, 2, 0, n_sigma);
end% if



%
% special for spectrometer measurements
%
if(plot_vector(8) == 1)
  n_plot = n_plot + 1;
  subplot(2,N_plot_cols, n_plot);
  % special case: CHER SPETRO ANA MAP Y to E
  load ~/Dropbox/temp/yE_nom.mat;
  pp(:,8) = my_2d_mapping(E_nom, y_nom, pp(:,2)); % E meas
%  pp(:,8)
  pp(:,9) = (pp(:,8) - pp(:,6)) ./ pp(:,6) * 100; % rel E error
  % n_perc
  n_perc = 10; % [%]
  p_err_ab_nperc = sum(pp(:,9) > n_perc) / length(pp(:,9)) * 100 
  my_densplot(pp, 6, 9, 0, n_sigma);
  subplot(2,N_plot_cols, n_plot+N_plot_cols);
  my_densplot(pp, 1, 8, 0, n_sigma); 
end% if


if(sum(plot_vector) > 0)
subplot(2,N_plot_cols,2);
  title([ '\mu_E [GeV] =' num2str(muE, '%.2f') ','...
          '\sigma_E / E [-] =' num2str(sigEE, '%.2E') ','...
          '\sigma_x [um]=' num2str(sigx, '%.1f') ','...
          '\sigma_y [um]=' num2str(sigy, '%.1f') ','...
          '\sigma_z [um]=' num2str(sigz, '%.1f') ','...
          '\beta_{x}=' num2str(betax, '%.3f' ) ','...
          '\alpha_{x}=' num2str(alphax, '%.3f' ) ','...
          '\beta_{y}=' num2str(betay, '%.3f' ) ','...
          '\alpha_{y}=' num2str(alphay, '%.3f' ) ','...
          '\epsilon_{Nx} [um] =' num2str(emnx, '%.2f') ','...
          '\epsilon_{Ny} [um] =' num2str(emny, '%.2f') ','...         
%          'c_{zp}=' num2str(corzp, '%.1E') ','...
%          'c_{zx}=' num2str(corzx, '%.1E') ...
        ]);
end% if

if(sum(plot_vector) > 0)
  set(0,'defaultaxesfontsize',18);
end% if
  

%
% beam slice quantities
%

% sort beam first
[Y, I] = sort(pp(:,3));
pp = pp(I, :);
[slice.mean_x, slice.sigma_x, slice.z, slice.N_z] = my_get_slice_var(pp, 3, 1);
[slice.mean_xp, slice.sigma_x, slice.z, slice.N_z] = my_get_slice_var(pp, 3, 4);
[slice.mean_y, slice.sigma_y, slice.z, slice.N_z] = my_get_slice_var(pp, 3, 2);
[slice.mean_E, slice.sigma_E, slice.z, slice.N_z] = my_get_slice_var(pp, 3, 6);
[slice.emnx,slice.emny,slice.betax,slice.alphax,slice.betay,slice.alphay] = my_get_slice_twiss(pp, slice.z);
0
if(plot_vector(7))  
  my_ana_beam_slice( slice );
end% if
  