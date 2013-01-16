%
% look at sliced beam quant's
%
if(1)

%load /Users/eadli/Dropbox/SLAC/quickpic/myData/n17_rampout_PI.mat
  
% show beam slices
n_beam = 1;
do_this_plot = 1;
for(n_3D_counter=1:length(qp)-1)
  s_timestep(n_3D_counter) = qp(n_3D_counter).s_timestep;
  slice.z = qp(n_3D_counter).PP(n_beam).slice.z;
  slice.N_z = qp(n_3D_counter).PP(n_beam).slice.N_z;
  n_non_zero = min(find(slice.N_z > 0))+1;
  slice.N_z = slice.N_z(n_non_zero:end);
  sigma_x(n_3D_counter) = qp(n_3D_counter).PP(n_beam).sigma_x;
  slice.mean_x = qp(n_3D_counter).PP(n_beam).slice.mean_x(n_non_zero:end);
  slice.sigma_x = qp(n_3D_counter).PP(n_beam).slice.sigma_x(n_non_zero:end);
  slice.mean_E = qp(n_3D_counter).PP(n_beam).slice.mean_E(n_non_zero:end);
  slice.sigma_E = qp(n_3D_counter).PP(n_beam).slice.sigma_E(n_non_zero:end);
  slice.betax = qp(n_3D_counter).PP(n_beam).slice.betax(n_non_zero:end);
  slice.betay = qp(n_3D_counter).PP(n_beam).slice.betay(n_non_zero:end);
  slice.emnx = qp(n_3D_counter).PP(n_beam).slice.emnx(n_non_zero:end);
  slice.emny = qp(n_3D_counter).PP(n_beam).slice.emny(n_non_zero:end);
  %sigma_x(n_3D_counter) = qp(n_3D_counter).PP(n_beam).sigma_x;

  if( do_this_plot)
  subplot(2,2,1);
  n_take_out_avg = 5:8;
  plot(slice.z(1:length(slice.mean_x)), slice.mean_x - mean(slice.mean_x(n_take_out_avg)), '-x');
  %xlabel('z [um]');
  ylabel('<x> [um]');
  grid on;
  title(['s=' num2str(s_timestep(n_3D_counter)*100, 3) ' [cm]']);
  %axis([-250 200 -10 10]);
  %
  subplot(2,2,2);
  plot(slice.z(1:length(slice.mean_x)), slice.sigma_x, '-x');
  %xlabel('z [um]');
  ylabel('\sigma_x [um]');
  grid on;
  %  
  subplot(2,2,2);
  plot(slice.z(1:length(slice.mean_x)), slice.mean_E, '-x');
  %xlabel('z [um]');
  ylabel('<E> [GeV]');
  grid on;
  %  
  subplot(2,2,1);
  plot(slice.z(1:length(slice.mean_x)), slice.emnx, '-x');
  %xlabel('z [um]');
  ylabel('\epsilon_{N,x} [um]');
  grid on;
  %  
  subplot(2,2,4);
  bar(slice.z(1:length(slice.mean_x)), slice.N_z);
  xlabel('z [um]');
  ylabel('\lambda_z [arb. units]');
  grid on;
  %  
  %subplot(2,2,4);
  %plot(slice.mean_x, slice.mean_E, '-x');
  %xlabel('<x> [um]');
  %ylabel('<E> [GeV]');
  %grid on;
  %
  % betax
  subplot(2,2,3);
  plot(slice.z(1:length(slice.mean_x)), slice.betax, '-x');
  xlabel('z [um]');
  ylabel('\beta_x [m]');
  grid on;
 
  subplot(2,2,1);
  title(['s=' num2str(s_timestep(n_3D_counter)*100, 3) ' [cm]']);
    pause;
  
  end% if

  %slice.mean_x_end(n_3D_counter) = mean(slice.mean_x(end-2:end-0));
  %slice.sigma_x_mean(n_3D_counter) = mean(slice.mean_x(end-2:end-0));
  n_end = length(slice.mean_x);
  n_start = 2;
%  n_mid = round(length(slice.mean_x)/2);
  n_mid = length(slice.mean_x)-5;
  slice.mean_x_end(n_3D_counter) = slice.mean_x(n_end)- 1*mean(slice.mean_x(n_take_out_avg));
  slice.mean_E_end(n_3D_counter) = slice.mean_E(n_end);
  %slice.mean_x_end(n_3D_counter) = mean(slice.mean_x(end-5));
  slice.sigma_x_mean(n_3D_counter) = mean(slice.sigma_x(n_start));
  slice.mean_E_mean(n_3D_counter) = mean(slice.mean_E(8));
  slice.sigma_E_mean(n_3D_counter) = mean(slice.sigma_E(8));
  %
  % compare with mid pos
  slice.mean_x_mid(n_3D_counter) = slice.mean_x(n_mid)- 1*mean(slice.mean_x(n_take_out_avg));
  slice.mean_E_mid(n_3D_counter) = slice.mean_E(n_mid);
  slice.mean_x_start(n_3D_counter) = slice.mean_x(n_start);
  slice.mean_E_start(n_3D_counter) = slice.mean_E(n_start);
  % twiss along
  slice.betax_mean(n_3D_counter) = mean(slice.betax(15:end));
  slice.betay_mean(n_3D_counter) = mean(slice.betay(15:end));
end% for

% use this to calc' fraction of charge above a certain energy
sum(slice.N_z(end-2:end))  / sum (slice.N_z)

%
% plot transverse parameters (hosing, erosion)
%
  subplot(2,1,1);
%  semilogy(s_timestep*1e2, abs(slice.mean_x_end), '-xb');
 semilogy(s_timestep*1e2, abs(slice.mean_x_end - 0*my_smooth(slice.mean_x_end, 100)))
%  hold on;
%  semilogy(s_timestep*1e2, abs(slice.mean_x_mid), '-xr');
  hold off;
legend(['<x> at z[um]=' num2str(slice.z(n_end), '%.0f') ',E_{final}[GeV]=' num2str(slice.mean_E_end(end),'%.0f')], ['<x> at z[um]=' num2str(slice.z(n_mid), '%.0f') ',E_{final}[GeV]=' num2str(slice.mean_E_mid(end),'%.0f')]);
  ylabel('<x> [um]'); % _{tail}
  xlabel('s [cm]');
axis([0 100 1e-2 1e2]);
grid on;
  %
    subplot(2,1,2);
  plot(s_timestep*1e2, slice.sigma_x_mean, '-x');
  hold on;
  plot(s_timestep*1e2, sigma_x, '-r');
  hold off;
  legend(['slice \sigma_x at z[um]=' num2str(slice.z(n_start), '%.0f') ',E_{final}[GeV]=' num2str(slice.mean_E_end(n_start),'%.0f')], 'bunch \sigma_x');
  ylabel('\sigma_x [um]'); % _{,head}
  xlabel('s [cm]');
  grid on;
%  grid on;
  %title(['s=' num2str(s_timestep*100, 3) ' [cm].   Time step: ' num2str(n_3D_timestep) ' [DT]']);

%
% plot energy gain for head and tail
%
  plot(s_timestep*1e2, abs(slice.mean_E_mean), '-x');
  ylabel('<E> [GeV]');
  xlabel('s [cm]');
  grid on;
  %


  %
% plot beta evolution
%
  delta_s = 0.01; % s in sim from which s started
  clf;
  plot(s_timestep*1e2-delta_s*1e2, slice.betax_mean, '-xb');
  hold on;
  plot(s_timestep*1e2-delta_s*1e2, slice.betay_mean, '-or');
  hold off;
  ylabel('\beta [m]');
  xlabel('s [cm]');
  grid on;
  %


%
% here, data for real plasma ramp to compare
%
  load /Users/eadli/Dropbox/SLAC/quickpic/myData/ideal_ramp_data_n1e17_z_beta.mat
  hold on;
  plot((z(end-63:end)-0.1+0.0016*2)*100, flipud(beta_evol(1:64)), '-xk');
  hold off;

  myaxis = axis;
  axis([0 10 myaxis(3) myaxis(4)])
  legend('\beta_{x,qp}', '\beta_{y,qp}', '\beta_{ideal}');
  
end% if




%
%
%  start CTM sliced analyzis
%
%

load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3_tilt0.mat;
load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3_tilt1.mat;
qp = qp(1:270); % cut when sim is eroded (45 cm)
%load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3_tilt5.mat;
load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3_tilt5_xp.mat;
qp = qp(1:238); % cut when sim hits box
load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3_tilt5_ExB.mat;
load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3_tilt50.mat;
qp = qp(1:101); % cut when sim hits box

%load -mat ~/Dropbox/SLAC/quickpic/myData/n3_tilt5_PI.mat 
load -mat ~/Dropbox/SLAC/quickpic/myData/n3_tilt5_se_PI.mat 

load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3_tilt0.mat;

% don't show slices with very small fraction of charge (outlier behaviour)
n_frac_acceptance = 1.0; % req at least this percentage of slice
n_good = (qp(1).PP(1).slice.N_z / sum(qp(1).PP(1).slice.N_z) * 100) > n_frac_acceptance;
myrange = min(find(n_good)) : max(find(n_good))
%n_z_range = min(find(n_good)) : max(find(n_good)) -1

% pre-scale all axes to extremes
myaxis.z_min = +1e99;
myaxis.z_max = -1e99;
myaxis.mean_x_min = +1e99;
myaxis.mean_x_max = -1e99;
myaxis.mean_xp_min = +1e99;
myaxis.mean_xp_max = -1e99;
myaxis.alphax_min = +1e99;
myaxis.alphax_max = -1e99;
myaxis.betax_min = +1e99;
myaxis.betax_max = -1e99;
myaxis.emnx_min = +1e99;
myaxis.emnx_max = -1e99;
myaxis.mean_E_min = +1e99;
myaxis.mean_E_max = -1e99;
myaxis.sigma_x_min = +1e99;
myaxis.sigma_x_max = -1e99;
for ix=1:size(qp,2),
  myaxis.z_min = min(myaxis.z_min, min(qp(ix).PP(1).slice.z(myrange)));
  myaxis.z_max = max(myaxis.z_max, max(qp(ix).PP(1).slice.z(myrange)));
  myaxis.mean_x_min = min(myaxis.mean_x_min, min(qp(ix).PP(1).slice.mean_x(myrange)));
  myaxis.mean_x_max = max(myaxis.mean_x_max, max(qp(ix).PP(1).slice.mean_x(myrange)));
  myaxis.mean_xp_min = min(myaxis.mean_xp_min, min(qp(ix).PP(1).slice.mean_xp(myrange)));
  myaxis.mean_xp_max = max(myaxis.mean_xp_max, max(qp(ix).PP(1).slice.mean_xp(myrange)));
  myaxis.alphax_min = min(myaxis.alphax_min, min(qp(ix).PP(1).slice.alphax(myrange)));
  myaxis.alphax_max = max(myaxis.alphax_max, max(qp(ix).PP(1).slice.alphax(myrange)));
  myaxis.betax_min = min(myaxis.betax_min, min(qp(ix).PP(1).slice.betax(myrange)));
  myaxis.betax_max = max(myaxis.betax_max, max(qp(ix).PP(1).slice.betax(myrange)));
  myaxis.emnx_min = min(myaxis.emnx_min, min(qp(ix).PP(1).slice.emnx(myrange)));
  myaxis.emnx_max = max(myaxis.emnx_max, max(qp(ix).PP(1).slice.emnx(myrange)));
  myaxis.mean_E_min = min(myaxis.mean_E_min, min(qp(ix).PP(1).slice.mean_E(myrange)));
  myaxis.mean_E_max = max(myaxis.mean_E_max, max(qp(ix).PP(1).slice.mean_E(myrange)));
  myaxis.sigma_x_min = min(myaxis.sigma_x_min, min(qp(ix).PP(1).slice.sigma_x(myrange)));
  myaxis.sigma_x_max = max(myaxis.sigma_x_max, max(qp(ix).PP(1).slice.sigma_x(myrange)));
end % for

for ix=1:size(qp,2),
% my_ana_beam_slice( qp(ix).PP(1).slice, myaxis, myrange, qp(ix).s_timestep );
%  pause(0.001);
%  pause;
  %n_mean_head = myrange(1);  % PI look at beam quantities in this part of the bunch
  %n_mean_tail = myrange(8);  % PI look at beam quantities in this part of the bunch
  n_mean_head = myrange(1);  % FI look at beam quantities in this part of the bunch
  n_mean_tail = myrange(11);  % FI look at beam quantities in this part of the bunch
  n_mean_tail = myrange(16);  % HOSING look at beam quantities in this part of the bunch
  s(ix) = qp(ix).s_timestep;
  s_mean_x_head(ix) = mean(qp(ix).PP(1).slice.mean_x(n_mean_head));
  s_mean_x_tail(ix) = mean(qp(ix).PP(1).slice.mean_x(n_mean_tail));
  s_mean_E_tail(ix) = mean(qp(ix).PP(1).slice.mean_E(n_mean_tail));
  s_mean_xp_tail(ix) = mean(qp(ix).PP(1).slice.mean_xp(n_mean_tail));
  s_tau_ion_n(ix) = min(find(qp(ix).PP(1).slice.sigma_x<20));
  s_tau_ion_z(ix) = qp(ix).PP(1).slice.z(s_tau_ion_n(ix));
  %s_ExB(ix) = qp(ix).PP(1).centr_cell_mag_v; % only for data with ExB activated
  s_mean_E_head(ix) = mean(qp(ix).PP(1).slice.mean_E(n_mean_head));
  s_mean_E_tail(ix) = mean(qp(ix).PP(1).slice.mean_E(n_mean_tail));
  omega_p = sqrt(qp(1).n0*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
  lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
  k_p = 2*pi/lambda_p;
  k_beta_tail = k_p / sqrt(2*s_mean_E_tail(ix)/0.511e-3);
  l_beta_tail(ix) = 2 * pi / k_beta_tail;
end% for
s_tau_ion_z = s_tau_ion_z - s_tau_ion_z(1);

%
% plotting
%

s_5 = s;
s_mean_x_head_5 = s_mean_x_head;
s_mean_x_tail_5 = s_mean_x_tail;


%load -mat /Users/eadli/Dropbox/SLAC/quickpic/myData/temp_tau_ion_conv.mat;
%plot(s*100, -conv_tau_ion*20, '-xk');
%plot(s_FI*100, s_tau_ion_z_FI, '-xr');
%hold on;
%plot(s*100, s_tau_ion_z_PI, '-xk');
%hold off;

set(0,'defaultaxesfontsize',28);

c_zx = 0.01;
plot(s_1*100, my_smooth(s_mean_x_head_1(1:length(s_1)), 100) / c_zx, '-.m');
hold on;
plot(s_1*100, s_mean_x_tail_1(1:length(s_1)) / c_zx, '-om');
c_zx = 0.05;
plot(s_5*100, my_smooth(s_mean_x_head_5(1:length(s_5)), 100) / c_zx, '-.b');
plot(s_5*100, s_mean_x_tail_5(1:length(s_5)) / c_zx, '-xb');
c_zx = 0.50;
plot(s_46*100, my_smooth(s_mean_x_head_46(1:length(s_46)), 100) / c_zx, '-.r');
plot(s_46*100, s_mean_x_tail_46(1:length(s_46)) / c_zx, '-+r');
hold off;
grid on;
xlabel('s [cm]');
ylabel('<xp> / c_{zx} [mrad]');
ylabel('<x> / c_{zx} [um]');
%ylabel('\tau_{ion} [um]');
grid on;
legend('c_{zx} = 0.01 rad,  <x>_{head}', ' c_{zx} = 0.01 rad, <x>_{tail}', 'c_{zx} = 0.05 rad,  <x>_{head}', ' c_{zx} = 0.05 rad, <x>_{tail}', 'c_{zx} = 0.46 rad,  <x>_{head}', ' c_{zx} = 0.46 rad, <x>_{tail}');
%legend('c_{zx} = 0.05 FI,  <x>_{head}', ' c_{zx} = 0.05 FI, <x>_{tail}', 'c_{zx} = 0.05 PI,  <x>_{head}', ' c_{zx} = 0.05 PI, <x>_{tail}');
%legend('c_{zx} = 0.05, <x>_{tail}', '\propto \int ds s x \tau_{sim}')
%legend('FI', 'PI')
myaxis = axis;
axis([0 40 -2500 500]);
title('CTM at head and tail, field-ionized plasma');
%
% XP
%
%plot(s*100, s_mean_xp_tail /1e3 / c_zx, '-xr');

%
% ENERGY
%
plot(s*100, s_mean_E_head, '-.xr');
hold on;
plot(s*100, s_mean_E_tail, '-xr');
hold off;
grid on;
xlabel('s [cm]');
ylabel('<E> [GeV]');
legend('<E>_{head}', '<E>_{tail}');
%
% BETA
%
plot(s*100, l_beta_tail*1e2, '-.xr');
grid on;
xlabel('s [cm]');
ylabel('\lambda_{\beta},{tail} [cm]');
%legend('<E>_{head}', '<E>_{tail}');


%
%
%

plot(s*100, s_ExB / c_zx, '-xr');
grid on;
hold on;
hold off;
xlabel('s [cm]');
ylabel('E x B / |B|^2 / c_{zx} [m/s]');
title('ExB drift at beam centroid / c_{zx}');


% summary for experiment
dDx_dczx = 2 / 0.4; % [mm / %]
x_s30_dczx = 1.1 / 0.05; % [mm / %]
xp_s30_dczx = 6.2 / 0.05; % mrad / %]

c_zx_range = 0:0.001:0.5;
Dx_range = 0:0.01:1;

subplot(2,1,1); 
plot(c_zx_range, x_s30_dczx * c_zx_range);
xlabel('c_{zx} [-]    (tilt)');
%plot(Dx_range, x_s30_dczx * Dx_range / dDx_dczx);
%xlabel('Dx [mm]    (dispersion at IP)');
ylabel('<x>   at plasma exit [mm]');
grid on;
subplot(2,1,2); 
plot(c_zx_range, xp_s30_dczx * c_zx_range);
xlabel('c_{zx} [-]    (tilt)');
%plot(Dx_range, xp_s30_dczx * Dx_range / dDx_dczx);
%xlabel('Dx [mm]    (dispersion at IP)');
ylabel('<xp>   at plasma exit [mrad]');
grid on;





%
%
%  start NJP emittance growth analysis
%
%

load ~/Dropbox/SLAC/quickpic/myData/NJP_0um.mat
qp_0 = qp;
load ~/Dropbox/SLAC/quickpic/myData/NJP_1um.mat
qp_1 = qp;
load ~/Dropbox/SLAC/quickpic/myData/NJP_10um.mat
qp_10 = qp;


plot(qp(end).PP(2).slice.z(1:end-1), qp_0(end).PP(2).slice.emnx, '-xk')
hold on;
plot(qp(end).PP(2).slice.z(1:end-1), qp_1(end).PP(2).slice.emnx, '-xr')
plot(qp(end).PP(2).slice.z(1:end-1), qp_10(end).PP(2).slice.emnx, '-xb')
hold off;


my_ana_beam_slice(qp_10(end).PP(2).slice);
