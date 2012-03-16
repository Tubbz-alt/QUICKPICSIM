%  function specially designed for quickpicsim_anabeam.m
function [gauss_sigma_z, gauss_mu_z, gauss_A_z] = my_gaussplot(fit_X, fit_Y, x_scaled, do_plot)

if nargin<4, do_plot = 1; end

fit_X = fit_X/1e6; % above scaled to 1e6 to get [um] om abscissa..
% gauss fit [from hist]
save -mat /tmp/fitfunc.dat fit_X fit_Y
N_guess_dividor = 1e5; % ~ # of particles ....
init_guess = [1e4/N_guess_dividor mean(fit_X) std(fit_X)];
result = fminsearch('mygaussfit2', init_guess);
gauss_sigma_z = abs(result(3));
gauss_mu_z = result(2);
gauss_A_z = result(1); % /gauss_sigma_z/sqrt(2*pi)
f_x = gauss_A_z*exp(-(0.5).*((fit_X-gauss_mu_z)/gauss_sigma_z).^2)';
if(do_plot)
  title_text = ['\sigma_{gauss} [um]=' num2str(gauss_sigma_z*1e6, '%.1f')];
  if(x_scaled)
    title_text(findstr(title_text, 'um') : findstr(title_text, 'um')+1) = 'mm';  
  end% if
  hold on;
  plot(fit_X*1e6, f_x, '-r');
  hold off;
  title(title_text);
end% if

