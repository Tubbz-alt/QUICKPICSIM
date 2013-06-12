clear all;
clf;

% ramp length
%z_ramp = 0.184; % [m]
z_ramp = 0.10; % [m]

n = 0;
p_min = 10;
p_max = 19;
n0_range = logspace(p_min, p_max, 100);

% emittance?  irrelevant

%gamma = 39139
gamma = 48924
gamma = 1e3/.511
gamma = 20e3/.511

% core plasma density
n0_range = 6.0e16; % osiris
n0_range = 2.211e16;
n0_range = 1.0e17;
n0_range = 2.5e17;

n0_range = 1e17;

n0_range = 1.0e16:1.0e16:5.0e17;

n0_range = 1e17;
gamma = 39139;
%n0_range = 2e16;
%gamma = 25e3/.511;

% intial beta error
d_beta_beta0 = 0.0;

% roll
roll_up = 0; % 1: up ramp, 0: down ramp

do_plot = 1;
for n0 = n0_range,
  n = n + 1;
  [b_mat(n), a_mat(n), b_0(n), a_0(n), M_tot, b_min(n), a_min(n), b_max(n), a_max(n), b_beat_rel(n)] = my_calc_ramp(n0, gamma, d_beta_beta0, do_plot, z_ramp, roll_up);
end% for

% vacuum waist
b_star = b_0 ./ (1+a_0.^2)
ds = a_0 .* b_star

% ramp demag
b_demag = b_star ./ b_mat

stop

semilogx(n0_range, b_demag, '-ok');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\beta^* / \beta_{mat}');

semilogx(n0_range, b_star*100, '-ok');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\beta^* [cm]');

semilogx(n0_range, ds*100, '-ok');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\Delta s_{ramp} [cm]');


stop

%
% 3 params graph
%
p_min = 15;
p_max = 18;
n0_range = logspace(p_min, p_max, 30);
d_beta_beta0 = 0;

gamma = 20e3/.511
z_ramp = 0.1; % [m]
do_plot = 0;
n = 0;
for n0 = n0_range,
  n = n + 1;
  [b_mat(n), a_mat(n), b_0(n), a_0(n), b_min(n), a_min(n), b_max(n), a_max(n), b_beat_rel(n)] = my_calc_ramp(n0, gamma, d_beta_beta0, do_plot, z_ramp);
end% for
% vacuum waist
b_star = b_0 ./ (1+a_0.^2);
ds = a_0 .* b_star;
% ramp demag
b_demag = b_star ./ b_mat;


subplot(1,2,2);
semilogx(n0_range, b_demag, '-ok');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\beta^* / \beta_{mat}');
hold on;
subplot(1,2,1);
semilogx(n0_range, b_star*100, '-ok');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\beta^* [cm]');
hold on;

% double length
gamma = 20e3/.511
z_ramp = 0.2; % [m]
do_plot = 0;
n = 0;
for n0 = n0_range,
  n = n + 1;
  [b_mat(n), a_mat(n), b_0(n), a_0(n), b_min(n), a_min(n), b_max(n), a_max(n), b_beat_rel(n)] = my_calc_ramp(n0, gamma, d_beta_beta0, do_plot, z_ramp);
end% for
% vacuum waist
b_star = b_0 ./ (1+a_0.^2);
ds = a_0 .* b_star;
% ramp demag
b_demag = b_star ./ b_mat;
subplot(1,2,2);
semilogx(n0_range, b_demag, '-+r');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\beta^* / \beta_{mat}');
subplot(1,2,1);
semilogx(n0_range, b_star*100, '-+r');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\beta^* [cm]');

% double energy
gamma = 40e3/.511
z_ramp = 0.1; % [m]
do_plot = 0;
n = 0;
for n0 = n0_range,
  n = n + 1;
  [b_mat(n), a_mat(n), b_0(n), a_0(n), b_min(n), a_min(n), b_max(n), a_max(n), b_beat_rel(n)] = my_calc_ramp(n0, gamma, d_beta_beta0, do_plot, z_ramp);
end% for
% vacuum waist
b_star = b_0 ./ (1+a_0.^2);
ds = a_0 .* b_star;
% ramp demag
b_demag = b_star ./ b_mat;
subplot(1,2,2);
semilogx(n0_range, b_demag, '-xb');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\beta^* / \beta_{mat}');
subplot(1,2,1);
semilogx(n0_range, b_star*100, '-xb');
grid on;
xlabel('n_0 [/cm^3]');
ylabel('\beta^* [cm]');

legend('E_0 = 20 GeV, z_{ramp} = 10 cm', 'E_0 = 20 GeV, z_{ramp} = 20 cm', 'E_0 = 40 GeV, z_{ramp} = 10 cm');




stop

%loglog(n0_range, b_mat, '-ok');
semilogx(n0_range, b_demag, '-ok');
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
