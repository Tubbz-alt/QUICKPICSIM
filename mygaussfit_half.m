function S = mygaussfit_half(params)

load -mat ~/Dropbox/temp/fitfunc.dat

x = fit_X;
y = fit_Y;

% Gaussian
%   f_x = A*exp(-0.5*(x/s)^2)
A = params(1);
mu = params(2);
s = params(3);
%f_x = A*(1/(s.*sqrt(2*pi)))*exp(-(0.5).*((x-mu)/s).^2);
f_x = A*exp(-(0.5).*((x-mu)/s).^2);
S = sum( ( f_x - y ).^2 );

% visu (comment when running)
if(1)
  axis([min(x) max(x) min(y) max(y)]);
  plot(x, f_x, '-xr');
  hold on;
  plot(x, y, '-xb');
  hold off;
%  legend('gauss fit', 'dist');
end% if

