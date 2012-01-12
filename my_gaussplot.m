%  function specially designed for quickpicsim_anabeam.m
function [gauss_sigma_z, gauss_mu_z, gauss_A_z] = my_gaussplot(fit_X, fit_Y, do_plot)

if nargin<3, do_plot = 1; end

fit_X = fit_X/1e6; % above scaled to 1e6 to get [um] om abscissa..
% gauss fit [from hist]
save -mat /tmp/fitfunc.dat fit_X fit_Y
init_guess = [1e2 0 std(fit_X)];
result = fminsearch('mygaussfit2', init_guess);
gauss_sigma_z = abs(result(3));
gauss_mu_z = result(2);
gauss_A_z = result(1); % /gauss_sigma_z/sqrt(2*pi)
f_x = gauss_A_z*exp(-(0.5).*((fit_X-gauss_mu_z)/gauss_sigma_z).^2)';
if(do_plot)
  hold on;
  plot(fit_X*1e6, f_x, '-r');
  hold off;
end% if