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
  slice_z = qp(n_3D_counter).PP(n_beam).slice_z;
  slice_N_z = qp(n_3D_counter).PP(n_beam).slice_N_z;
  n_non_zero = min(find(slice_N_z > 0))+1;
  slice_N_z = slice_N_z(n_non_zero:end);
  sigma_x(n_3D_counter) = qp(n_3D_counter).PP(n_beam).sigma_x;
  slice_mean_x = qp(n_3D_counter).PP(n_beam).slice_mean_x(n_non_zero:end);
  slice_sigma_x = qp(n_3D_counter).PP(n_beam).slice_sigma_x(n_non_zero:end);
  slice_mean_E = qp(n_3D_counter).PP(n_beam).slice_mean_E(n_non_zero:end);
  slice_sigma_E = qp(n_3D_counter).PP(n_beam).slice_sigma_E(n_non_zero:end);
  slice_betax = qp(n_3D_counter).PP(n_beam).slice_betax(n_non_zero:end);
  slice_betay = qp(n_3D_counter).PP(n_beam).slice_betay(n_non_zero:end);
  slice_emnx = qp(n_3D_counter).PP(n_beam).slice_emnx(n_non_zero:end);
  slice_emny = qp(n_3D_counter).PP(n_beam).slice_emny(n_non_zero:end);
  %sigma_x(n_3D_counter) = qp(n_3D_counter).PP(n_beam).sigma_x;

  if( do_this_plot)
  subplot(2,2,1);
  n_take_out_avg = 5:8;
  plot(slice_z(1:length(slice_mean_x)), slice_mean_x - mean(slice_mean_x(n_take_out_avg)), '-x');
  %xlabel('z [um]');
  ylabel('<x> [um]');
  grid on;
  title(['s=' num2str(s_timestep(n_3D_counter)*100, 3) ' [cm]']);
  %axis([-250 200 -10 10]);
  %
  subplot(2,2,2);
  plot(slice_z(1:length(slice_mean_x)), slice_sigma_x, '-x');
  %xlabel('z [um]');
  ylabel('\sigma_x [um]');
  grid on;
  %  
  subplot(2,2,2);
  plot(slice_z(1:length(slice_mean_x)), slice_mean_E, '-x');
  %xlabel('z [um]');
  ylabel('<E> [GeV]');
  grid on;
  %  
  subplot(2,2,1);
  plot(slice_z(1:length(slice_mean_x)), slice_emnx, '-x');
  %xlabel('z [um]');
  ylabel('\epsilon_{N,x} [um]');
  grid on;
  %  
  subplot(2,2,4);
  bar(slice_z(1:length(slice_mean_x)), slice_N_z);
  xlabel('z [um]');
  ylabel('\lambda_z [arb. units]');
  grid on;
  %  
  %subplot(2,2,4);
  %plot(slice_mean_x, slice_mean_E, '-x');
  %xlabel('<x> [um]');
  %ylabel('<E> [GeV]');
  %grid on;
  %
  % betax
  subplot(2,2,3);
  plot(slice_z(1:length(slice_mean_x)), slice_betax, '-x');
  xlabel('z [um]');
  ylabel('\beta_x [m]');
  grid on;

  subplot(2,2,1);
  title(['s=' num2str(s_timestep(n_3D_counter)*100, 3) ' [cm]']);
    pause;
  
  end% if

  %slice_mean_x_end(n_3D_counter) = mean(slice_mean_x(end-2:end-0));
  %slice_sigma_x_mean(n_3D_counter) = mean(slice_mean_x(end-2:end-0));
  n_end = length(slice_mean_x);
  n_start = 2;
%  n_mid = round(length(slice_mean_x)/2);
  n_mid = length(slice_mean_x)-5;
  slice_mean_x_end(n_3D_counter) = slice_mean_x(n_end)- 1*mean(slice_mean_x(n_take_out_avg));
  slice_mean_E_end(n_3D_counter) = slice_mean_E(n_end);
  %slice_mean_x_end(n_3D_counter) = mean(slice_mean_x(end-5));
  slice_sigma_x_mean(n_3D_counter) = mean(slice_sigma_x(n_start));
  slice_mean_E_mean(n_3D_counter) = mean(slice_mean_E(8));
  slice_sigma_E_mean(n_3D_counter) = mean(slice_sigma_E(8));
  %
  % compare with mid pos
  slice_mean_x_mid(n_3D_counter) = slice_mean_x(n_mid)- 1*mean(slice_mean_x(n_take_out_avg));
  slice_mean_E_mid(n_3D_counter) = slice_mean_E(n_mid);
  slice_mean_x_start(n_3D_counter) = slice_mean_x(n_start);
  slice_mean_E_start(n_3D_counter) = slice_mean_E(n_start);
  % twiss along
  slice_betax_mean(n_3D_counter) = mean(slice_betax(15:end));
  slice_betay_mean(n_3D_counter) = mean(slice_betay(15:end));
end% for

% use this to calc' fraction of charge above a certain energy
sum(slice_N_z(end-2:end))  / sum (slice_N_z)

%
% plot transverse parameters (hosing, erosion)
%
  subplot(2,1,1);
%  semilogy(s_timestep*1e2, abs(slice_mean_x_end), '-xb');
 semilogy(s_timestep*1e2, abs(slice_mean_x_end - 0*my_smooth(slice_mean_x_end, 100)))
%  hold on;
%  semilogy(s_timestep*1e2, abs(slice_mean_x_mid), '-xr');
  hold off;
legend(['<x> at z[um]=' num2str(slice_z(n_end), '%.0f') ',E_{final}[GeV]=' num2str(slice_mean_E_end(end),'%.0f')], ['<x> at z[um]=' num2str(slice_z(n_mid), '%.0f') ',E_{final}[GeV]=' num2str(slice_mean_E_mid(end),'%.0f')]);
  ylabel('<x> [um]'); % _{tail}
  xlabel('s [cm]');
axis([0 100 1e-2 1e2]);
grid on;
  %
    subplot(2,1,2);
  plot(s_timestep*1e2, slice_sigma_x_mean, '-x');
  hold on;
  plot(s_timestep*1e2, sigma_x, '-r');
  hold off;
  legend(['slice \sigma_x at z[um]=' num2str(slice_z(n_start), '%.0f') ',E_{final}[GeV]=' num2str(slice_mean_E_end(n_start),'%.0f')], 'bunch \sigma_x');
  ylabel('\sigma_x [um]'); % _{,head}
  xlabel('s [cm]');
  grid on;
%  grid on;
  %title(['s=' num2str(s_timestep*100, 3) ' [cm].   Time step: ' num2str(n_3D_timestep) ' [DT]']);

%
% plot energy gain for head and tail
%
  plot(s_timestep*1e2, abs(slice_mean_E_mean), '-x');
  ylabel('<E> [GeV]');
  xlabel('s [cm]');
  grid on;
  %


  %
% plot beta evolution
%
  delta_s = 0.01; % s in sim from which s started
  clf;
  plot(s_timestep*1e2-delta_s*1e2, slice_betax_mean, '-xb');
  hold on;
  plot(s_timestep*1e2-delta_s*1e2, slice_betay_mean, '-or');
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



