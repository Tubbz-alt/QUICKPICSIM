clear all;
clf;

% core plasma density
n0     = 2.3E13; % [cm^-3]
%n0     = 6.35E16 % [cm^-3]
gamma = 39139
% emittance?  irrelevant

n = 0;
p_min = 10;
p_max = 19;
n0_range = logspace(p_min, p_max, 100);

n0_range = 4.5e17;
n0_range = 2.3e16;
n0_range = 3.6e17;
n0_range = 5.4e16;
n0_range = 6.0e16; % osiris
n0_range = 1.0e17;
n0_range = 3.6e17;

% intial beta error
d_beta_beta0 = 0.0;

for n0 = n0_range,
  n = n + 1;
  [b_mat(n), a_mat(n), b_0(n), a_0(n), b_min(n), a_min(n), b_max(n), a_max(n), b_beat_rel(n)] = my_calc_ramp(n0, gamma, d_beta_beta0, 1);
end% for

stop

loglog(n0_range, b_mat, '-ok');
hold on;
%loglog(n0_range, b_max, '-g');
%loglog(n0_range, b_min, '-m');
loglog(n0_range, b_0, '-xm');
hold off;
xlabel('n_0 [cm^{-3}]');
ylabel('\beta [m]');

legend('\beta_{mat}', '\beta_{max}', '\beta_{min}', '\beta_{0}');
legend('\beta_{mat}', '\beta_{0}');

myaxis = axis;
axis([10^p_min  10^p_max  myaxis(3)  myaxis(4)]);

%grid on;

stop

%
% scale to an emittance of 30 um (predicated IP emittance for y)
%

emitt_N = 30e-6;
sigma_mat = sqrt(emitt_N / gamma * b_mat);
sigma_0 = sqrt(emitt_N / gamma * b_0);

loglog(n0_range, sigma_mat*1e6, '-ok');
hold on;
loglog(n0_range, sigma_0*1e6, '-xm');
hold off;
xlabel('n_0 [cm^{-3}]');
ylabel('\sigma [um]        (\epsilon_N = 30 um, \gamma = 39139)');
legend('\sigma_{mat}', '\sigma_{0}');
myaxis = axis;
axis([10^p_min  10^p_max  myaxis(3)  myaxis(4)]);
axis([10^p_min  10^p_max  1e-1 1e2]);
grid on;

stop

%
% plot relative beta beat
% 
semilogx(n0_range, b_beat_rel*1e2, '-ok');
xlabel('n_0 [cm^{-3}]');
ylabel('max. d\beta/\beta in plasma flat-top)  [%],   for d(\beta_0)/(\beta_0)=0.1');
myaxis = axis;
axis([10^p_min  10^p_max  myaxis(3)  myaxis(4)]);
grid on;
