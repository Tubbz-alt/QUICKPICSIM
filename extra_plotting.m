
%
% some extra plotting here
%

do_extra_plot = 0;

if(do_extra_plot)
 N_db = 3;
 N_wb = 1;
 
  mean_E = my_flatten_struct(qp, 'mean_E');
  std_E = my_flatten_struct(qp, 'std_E');
  abscissa = [qp(:).s_timestep];
  
  plot(abscissa, mean_E(:,1), '-xb');
  hold on;
  plot(abscissa, mean_E(:,2), '-or');
  plot(abscissa, mean_E(:,1)*N_db + mean_E(:,2)*N_wb, '-k'); % extract info here
  hold off;
  grid on;
  h_label = xlabel('s [cm]');
  set(h_label,'FontSize',14);
  h_label = ylabel('E_{part} [GeV],  E_{tot} [a.u.]');
  set(h_label,'FontSize',14);
  h_legend = legend('E_{mean,DB}', 'E_{mean,WB}', 'E_{tot,DB} + E_{tot,WB}');
  %set(h_legend,'FontSize',14);
  axis([0 0.8 0 100]);
% efficiency calcs  
 %E_final_db = qp(1).PP(1d).N_db;
 %E_final_wb = qp(1).PP(2).N_wb;
 E_final_db = mean_E(end,1)*N_db;
 E_final_wb = mean_E(end,2)*N_wb
 E_initial_db = mean_E(1,1)*N_db
 E_initial_wb = mean_E(1,2)*N_wb
 dE_db = E_final_db-E_initial_db
 dE_wb = E_final_wb-E_initial_wb
 final_trans_eff = dE_wb / dE_db
end% if

% convergense tests
if(0)
  load ~/temp/conv.mat;
  ZZ = [313/256:313/256:313];
  plot(ZZ, FEZ998_1(:,2^9/2+1), '-ob')
  hold on;
  %plot(ZZ, FEZ998_4(:,2^9/2+1), '-k')
  plot(ZZ, FEZ008_1(:,2^10/2+1), '-or')
  %plot(ZZ, FEZ008_2(:,2^10/2+1), '-k')
  plot(ZZ, FEZ118_1(:,2^11/2+1), '-og')
  %plot(ZZ, FEZ118_2(:,2^11/2+1), '-k')
  hold off;
  axis([0 313 -140 80]);
  grid on;
  xlabel('z [um]');
  ylabel('E_z [GV/m]');
  legend('INDX,Y,Z : 9x9x8', 'INDX,Y,Z : 10x10x8', 'INDX,Y,Z : 11x11x8');
  legend('9x9x8, NPx1', '9x9x8, NPx2');
  legend('10x10x8, NPx1', '10x10x8, NPx4');
  legend('11x11x8, NPx1', '11x11x8, NPx4');
end% 
  





% temp old stuff
  if(0)
% 

% plot total beam quantitites

s_timestep = [qp(1:end).s_timestep];
if(do_QEB_3D)
  Q_total = [qp(1:end).Q_total];
end% if
sigma_x = my_flatten_struct(qp, 'sigma_x');
mean_x = my_flatten_struct(qp, 'mean_x');
sigma_y = my_flatten_struct(qp, 'sigma_y');
sigma_z = my_flatten_struct(qp, 'sigma_z');
mean_y = my_flatten_struct(qp, 'mean_y');
emitt_x = my_flatten_struct(qp, 'emitt(1)');
emitt_y = my_flatten_struct(qp, 'emitt(2)');

plot(s_timestep*100, emitt_x*1e6);
xlabel('Timestep [cm]');
ylabel('\epsilon_N,x [um]');
grid on;
if( size(emitt_x, 2) == 2 );
  legend('DB', 'WB')
end; % if

plot(s_timestep*100, emitt_y*1e6);
xlabel('Timestep [cm]');
ylabel('\epsilon_N,y [um]');
grid on;
if( size(emitt_x, 2) == 2 );
  legend('DB', 'WB')
end; % if

plot(s_timestep*100, sigma_x);
xlabel('Timestep [cm]');
ylabel('\sigma_x [um]');
grid on;
if( size(emitt_x, 2) == 2 );
  legend('DB', 'WB')
end; % if


plot(s_timestep*100, sigma_y);
xlabel('Timestep [cm]');
ylabel('\sigma_y [um]');
grid on;
if( size(emitt_x, 2) == 2 );
  legend('DB', 'WB')
end; % if


plot(s_timestep*100, sigma_z);
xlabel('Timestep [cm]');
ylabel('\sigma_z [um]');
grid on;
if( size(emitt_x, 2) == 2 );
  legend('DB', 'WB')
end; % if

plot(s_timestep*100, mean_x);
xlabel('Timestep [cm]');
ylabel('<x> [um]');
grid on;
if( size(emitt_x, 2) == 2 );
  legend('DB', 'WB')
end; % if



plot(s_timestep*100, mean_y);
xlabel('Timestep [cm]');
ylabel('<y> [um]');
grid on;
if( size(emitt_x, 2) == 2 );
  legend('DB', 'WB')
end; % if



%plot(Q_total / max(Q_total), '-r')
%xlabel('Timestep [cm]');
%ylabel('Totat integrated charge QEB [a.u.]');
%grid on;


% either one beam or both beams
pp = qp(1).PP(1).BEAM;
%pp = qp_BEAMS;

subplot(2,2,1);
plot( pp(:,3)*1e0, pp(:,6)/1e0, 'o');
xlabel('z [um]'); 
ylabel('E [GeV]');
myaxis = axis;
axis([min(ZZ)  max(ZZ)  myaxis(3)  myaxis(4)]);
grid on;

subplot(2,2,3);
[NN, XX] = hist(pp(:,3)*1e0, 51);
bar(XX, NN);
xlabel('z [um]'); 
ylabel('count [a.u.]');
myaxis = axis;
axis([min(ZZ)  max(ZZ)  myaxis(3)  myaxis(4)]);
grid on;

subplot(2,2,2);
plot( pp(:,3)*1e0, pp(:,1)*1e0, 'o');
xlabel('z [um]'); 
ylabel('x [um]');
myaxis = axis;
axis([min(ZZ)  max(ZZ)  myaxis(3)  myaxis(4)]);
grid on;

subplot(2,2,4);
plot( pp(:,1)*1e0, pp(:,4)/1e0, 'o');
xlabel('x [um]'); 
ylabel('xp [urad]');
grid on;

%
%

[NN, XX] = hist(pp(:,1)*1e0, 102);
bar(XX, NN);
xlabel('x [um]'); 
ylabel('count [a.u.]');
grid on;

end% if 0

%
% look at sliced beam quant's
%
if(0)
  
% show beam slices
n_beam = 1;
do_this_plot = 1;
for(n_3D_counter=1:length(qp)-1)
  s_timestep(n_3D_counter) = qp(n_3D_counter).s_timestep;
  z_slices = qp(n_3D_counter).PP(n_beam).slice_z;
  slice_N_z = qp(n_3D_counter).PP(n_beam).slice_N_z;
  n_non_zero = min(find(slice_N_z > 0))+1;
  slice_N_z = slice_N_z(n_non_zero:end);
  sigma_x(n_3D_counter) = qp(n_3D_counter).PP(n_beam).sigma_x;
  slice_mean_x = qp(n_3D_counter).PP(n_beam).slice_mean_x(n_non_zero:end);
  slice_sigma_x = qp(n_3D_counter).PP(n_beam).slice_sigma_x(n_non_zero:end);
  slice_mean_E = qp(n_3D_counter).PP(n_beam).slice_mean_E(n_non_zero:end);
  slice_sigma_E = qp(n_3D_counter).PP(n_beam).slice_sigma_E(n_non_zero:end);
  %sigma_x(n_3D_counter) = qp(n_3D_counter).PP(n_beam).sigma_x;

  if( do_this_plot)
  subplot(2,2,1);
  n_take_out_avg = 5:8;
  plot(z_slices(1:length(slice_mean_x)), slice_mean_x - mean(slice_mean_x(n_take_out_avg)), '-x');
  %xlabel('z [um]');
  ylabel('<x> [um]');
  grid on;
  title(['s=' num2str(s_timestep(n_3D_counter)*100, 3) ' [cm]']);
  %axis([-250 200 -10 10]);
  %
  subplot(2,2,2);
  plot(z_slices(1:length(slice_mean_x)), slice_sigma_x, '-x');
  %xlabel('z [um]');
  ylabel('\sigma_x [um]');
  grid on;
  %  
  subplot(2,2,3);
  plot(z_slices(1:length(slice_mean_x)), slice_mean_E, '-x');
  xlabel('z [um]');
  ylabel('<E> [GeV]');
  grid on;
  %  
  subplot(2,2,4);
  bar(z_slices(1:length(slice_mean_x)), slice_N_z);
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
legend(['<x> at z[um]=' num2str(z_slices(n_end), '%.0f') ',E_{final}[GeV]=' num2str(slice_mean_E_end(end),'%.0f')], ['<x> at z[um]=' num2str(z_slices(n_mid), '%.0f') ',E_{final}[GeV]=' num2str(slice_mean_E_mid(end),'%.0f')]);
  ylabel('<x> [um]'); % _{tail}
  xlabel('s [cm]');
%axis([0 100 1e-2 1e2]);
grid on;
  %
    subplot(2,1,2);
  plot(s_timestep*1e2, slice_sigma_x_mean, '-x');
  hold on;
  plot(s_timestep*1e2, sigma_x, '-r');
  hold off;
  legend(['slice \sigma_x at z[um]=' num2str(z_slices(n_start), '%.0f') ',E_{final}[GeV]=' num2str(slice_mean_E_end(n_start),'%.0f')], 'bunch \sigma_x');
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

  
  

end% if

if(0)
%
% plot some quantitites along plasma
%
%[AX, H1, H2] = plotyy([qp.s_timestep]*100, [qp.FEZ_max], [qp.s_timestep]*100, [qp.QEB_max], 'plot');
%set(H1,'Color', 'r');
%set(H2,'Color', '-b');
%set(H1,'Color', 'r');
%set(H2,'Color', '-b');
load -mat ~/temp/match.mat
%load -mat ~/temp/4025match.mat
plot([qp.s_timestep]*100, [qp.FEZ_max] /  max([qp.FEZ_max]), '-xr');
hold on;
plot([qp.s_timestep]*100, [qp.QEB_max]/  max([qp.QEB_max]), '-xb');
%plot([qp.s_timestep]*100, [qp.QEB_max]/10, '-ok');
load -mat ~/temp/nomatch.mat
qp = qp(1:end-1);
max_field = max([qp.FEZ_max]);
plot([qp.s_timestep]*100, [qp.FEZ_max] /  max([qp.FEZ_max]), '-om');
plot([qp.s_timestep]*100, [qp.QEB_max]/  max([qp.QEB_max]), '-ok');
hold off;
grid on;
xlabel('s [cm]');
ylabel('E_z [GV/m]  and  n_b [10 n_{p0}]');
ylabel('E_z / E_{z,max}    ,    n_b / n_{b,max}');
legend('E_z, a) 20x20x20 um^3, hor. matched', 'n_{b}, a) 20x20x20 um^3, hor. matched', 'E_z, b) 20x20x20 um^3, no lens', 'n_{b}, b) 20x20x20 um^3, no lens');
%legend('E_z, a) 25x25x40 um^3, matched', 'n_{b}, a)  25x25x40 um^3, matched');
myaxis = axis;
axis([0 80 myaxis(3) myaxis(4)]);

end% if

%
% QP to ELE
%
% pp_qp = qp_BEAMS;
% pp = my_qp2ele(pp_qp);
% save -ascii -double ~/temp/dist_s50.asc pp
 