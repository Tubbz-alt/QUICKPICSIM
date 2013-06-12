z_ion_head =  my_flatten_struct(qp, 'z_ion_head');
x_ion_head =  my_flatten_struct(qp, 'x_ion_head');
x_zero_EcB =  my_flatten_struct(qp, 'x_zero_EcB');
EcB_zero =  my_flatten_struct(qp, 'EcB_zero');
EcB_at_beam =  my_flatten_struct(qp, 'EcB_at_beam');
plot([qp(:).s_timestep]*100, my_flatten_struct(qp, 'EcB_at_beam'), '-xk');
hold on;
plot([qp(:).s_timestep]*100, my_flatten_struct(qp, 'EcB_zero'), '-+r');
hold off;
grid on;
xlabel('s [cm]');
ylabel('E-cB [GV/m]');
legend('<EcB> beam',  '<EcB> grid with min. EcB', 'Location', 'Best');
title_text = datadir;
title(title_text, 'interpreter', 'None');


stop

plot([qp(:).s_timestep]*100, my_flatten_struct(qp, 'x_beam_EcB'), '-xk');
hold on;
plot([qp(:).s_timestep]*100, my_flatten_struct(qp, 'x_zero_EcB'), '-+r');
hold off;
grid on;
xlabel('s [cm]');
ylabel('<x> [um]');
legend('<x> Beam',  '<x> F\perp = 0', 'Location', 'Best');
title_text = datadir;
title(title_text, 'interpreter', 'None');

stop
%
%
%


% c_fudge
amp_corr = 5.0494; % crcy was wrong by a factor (k_p / sqrt(2))*z
PI_matched_hose = [0.10 0.27 0.63 0.95  1.15  1.40 1.51 1.63 1.64 1.64 1.64] / amp_corr
PI_unmatched_hose = [NaN 0.27 0.65  0.87  0.97  1.0  1.01 NaN  1.0 NaN 1.00] / amp_corr
FI_matched_hose = [NaN 0.25  0.39  0.64  0.86   0.94  1.06  1.12  1.03  1.03  1.01 ] / amp_corr;
FI_unmatched_hose =  [NaN 0.25  0.37  0.60   0.73  0.78  0.77 NaN  0.82 NaN 0.82] / amp_corr;
sr = 5./sqrt([2^99 2^18 2^14 2^12 2^8 2^5 2^4 2^3 2^2 2^1 2^0 2^(-1)]);
sr = sr(2:end);

plot(sr, fliplr(FI_matched_hose), '-xb');
hold on;
plot(sr, fliplr(FI_unmatched_hose), ':ob');
plot(sr, fliplr(PI_matched_hose), '-xr');
plot(sr, fliplr(PI_unmatched_hose), ':or');
hold off;
grid on;
xlabel('\sigma_r [um]');
ylabel('c_r \times c_\psi [-]');
legend('FI mat', 'FI unmat', 'PI mat', 'PI ummat', 'Location', 'Best');
 myaxis = axis;
 axis([1.23 myaxis(2) myaxis(3) myaxis(4)]);

%stop

%
% hosing growth analysis (corrected for ctm)
%
z_ion_head =  my_flatten_struct(qp, 'z_ion_head');
x_ion_head =  my_flatten_struct(qp, 'x_ion_head');
x_tail =  my_flatten_struct(qp, 'x_tail');
x_hose = x_tail-x_ion_head;
plot([qp(:).s_timestep]*100, x_hose, '-xb'); % osc
% abs osc 
ind_raise = 1;
ind_n = 2;
while ind_n <= length(x_hose),
  x_fitmax = 55; % before upset by sheath
  if( (abs(x_hose(ind_n)) > abs(x_hose(ind_raise(end)))) && (abs(x_hose(ind_n)) < x_fitmax) )
    ind_raise = [ind_raise ind_n];
  end% if
  if(abs(x_hose(ind_n)) > x_fitmax)
    ind_n = length(x_hose)+1;
  end% if
  ind_n = ind_n + 1;
end% while
%plot([qp(:).s_timestep]*100, abs(x_hose), '-xb'); % abs osc
hold on;
s_raise = [qp(ind_raise).s_timestep];
x_raise = abs(x_hose(ind_raise));
plot(s_raise*100, x_raise, '-or'); % abs osc
% normalize x_raise to x_raise(1) = 1;
%x_raise = x_raise - x_raise(1) + 1;
amp_sim = x_raise / x_raise(1);
% fit exp growth param
n0 = 1e17;
n_sigma_fit = (3+qp(1).PP(n_beam).hose_calc.n_sigma);
z = n_sigma_fit*qp(1).PP(1).sigma_z * 1e-6; % z for hose fit corresponds to z where tail is tracked
Gamma = qp_Gamma(1);
save -mat /tmp/fitfunc.dat n0 z s_raise amp_sim Gamma
addpath('~/Dropbox/SLAC/E200/E200_hose/');
fit_result = fminsearch('E200_hose_fit', 1.0);
c_fudge_fit = fit_result(1)
% plot fit
[A, amp_fit] = E200_calc_hose_A(n0, Gamma, s_raise, z, c_fudge_fit);
plot(s_raise*100, x_raise(1)*amp_fit, '-+k'); % abs osc
hold off
grid on;
xlabel('s [cm]');
ylabel('x_{tail} [um]');
legend(['x @ ' sprintf('%0.1f', n_sigma_fit) '\sigma_z'], 'envl.', ['envl. fit, c_{fudge} = ' sprintf('%0.1e', c_fudge_fit)], 'Location', 'Best');
title_text = datadir;
title(title_text, 'interpreter', 'None');


stop

%
% some extra plotting here
%
N_BEAM = 2; % which beam to look at (PWFALC)
%N_BEAM = 1; % which beam to look at (hosing)

% which slice to analyse
n_myslice = 20; 

mean_E = my_flatten_struct(qp, 'mean_E');
std_E = my_flatten_struct(qp, 'std_E');
mean_x = my_flatten_struct(qp, 'mean_x');
mean_xp = my_flatten_struct(qp, 'mean_xp');
mean_y = my_flatten_struct(qp, 'mean_y');
mean_yp = my_flatten_struct(qp, 'mean_yp');
emitt_x = my_flatten_struct(qp, 'emitt(1)');
emitt_y = my_flatten_struct(qp, 'emitt(2)');
abscissa = [qp(:).s_timestep];

plot(abscissa*100, mean_x(:,N_BEAM), '-xb');
hold on;
plot(abscissa*100, mean_xp(:,N_BEAM), '-or');
%plot(abscissa*100, mean_y(:,N_BEAM), '-xb');
%hold on;
%plot(abscissa*100, mean_yp(:,N_BEAM), '-or');
hold off;
xlabel('s [cm]');
ylabel('<x>, <xp> [um, urad]');
legend('<x> [um]', '<xp> [urad]');
%ylabel('<y>, <yp> [um, urad]');
%legend('<y> [um]', '<yp> [urad]');
grid on;

pause;

plot(abscissa*100, mean_E(:,N_BEAM), '-xb');
xlabel('s [cm]');
ylabel('<E> [GeV]');
grid on;

pause;

% look at last betax 
sliceinfo = my_flatten_struct(qp, 'slice');
N_3D = size(sliceinfo, 1);
betax = [sliceinfo(:, N_BEAM).betax];
betay = [sliceinfo(:, N_BEAM).betay];
emnx = [sliceinfo(:, N_BEAM).emnx];
emny = [sliceinfo(:, N_BEAM).emny];
sigma_x = [sliceinfo(:, N_BEAM).sigma_x];
sigma_y = [sliceinfo(:, N_BEAM).sigma_y];
if (do_ana_slice__gauss_fit)
  gaussx = [sliceinfo(:, N_BEAM).gaussx];
  gaussy = [sliceinfo(:, N_BEAM).gaussy];
end

N_slices = size(betax, 2) / N_3D;
n_myslice = round(N_slices/2)+1; % which slice to analyse
n_slice_fe = N_slices - n_myslice;
betax_slice = betax((N_slices-n_slice_fe):N_slices:(N_3D*N_slices-n_slice_fe));
betay_slice = betay((N_slices-n_slice_fe):N_slices:(N_3D*N_slices-n_slice_fe));
emnx_slice = emnx((N_slices-n_slice_fe):N_slices:(N_3D*N_slices-n_slice_fe));
emny_slice = emny((N_slices-n_slice_fe):N_slices:(N_3D*N_slices-n_slice_fe));
sigma_x_slice = sigma_x((N_slices-n_slice_fe):N_slices:(N_3D*N_slices-n_slice_fe));
%sigma_y_slice = sigma_y((N_slices-n_slice_fe):N_slices:(N_3D*N_slices-n_slice_fe));
if (do_ana_slice__gauss_fit)
  gaussx_slice = gaussx((N_slices-n_slice_fe):N_slices:(N_3D*N_slices-n_slice_fe));
  gaussy_slice = gaussy((N_slices-n_slice_fe):N_slices:(N_3D*N_slices-n_slice_fe));
end% if

plot(abscissa*100, betax_slice*100, '-x');
grid on;
axis([0 30 0.0 5]);
xlabel('s [cm]');
ylabel('\beta_x [cm]');



%
subplot(3,1,1);
plot(abscissa*1e2, emnx_slice, '-r');
hold on;
plot(abscissa*1e2, emny_slice, '-b');
hold off;
grid on;
xlabel('s [cm]');
ylabel('\epsilon_N [um]');
legend('x', 'y');
myaxis = axis;
%axis([0 40 myaxis(3) myaxis(4)]);
title(['Slice ' num2str(n_myslice) ' out of ' num2str(N_slices) '. z pos ' num2str(slice.z(N_BEAM, n_myslice), 3) ' um']);

subplot(3,1,2);
plot(abscissa*1e2, betax_slice, '-r');
hold on;
plot(abscissa*1e2, betay_slice, '-b');
hold off;
grid on;
xlabel('s [cm]');
ylabel('\beta [m]');
legend('x', 'y');

subplot(3,1,3);
plot(abscissa*1e2, sigma_x_slice, '-r');
hold on;
%plot(abscissa*1e2, sigma_y_slice, '-b');
hold off;
grid on;
xlabel('s [cm]');
ylabel('\sigma [um]');
legend('x', 'y');

stop
subplot(4,1,4);
plot(abscissa*1e2, gaussx_slice*1e6, '-r');
hold on;
plot(abscissa*1e2, gaussy_slice*1e6, '-b');
hold off;
grid on;
xlabel('s [cm]');
ylabel('\sigma_{gauss} [um]');
legend('x', 'y');



stop


%
% mean energy and plasma density
%
%s = s_ramp;
%n = n_ramp;
s = s_ft;
n = n_ft;
[AX,H1,H2] = plotyy(abscissa*100, mean_E(:,N_BEAM), s/1e4, n, 'plot');
xlabel('s [cm]');
grid on;
set(H1,'Color','b')
set(H2,'Color','g')
set(get(AX(1),'Ylabel'),'String','<E> [GeV]'); 
set(get(AX(2),'Ylabel'),'String', 'n_p / n_0 [-]');   
set(get(AX(1),'Ylabel'),'FontSize',22); 
set(get(AX(2),'Ylabel'),'FontSize',22); 
set(get(AX(1),'Ylabel'),'Color','blue'); 
set(get(AX(2),'Ylabel'),'Color','green'); 
%set(H1,'LineStyle','-')
%set(H2,'LineStyle','-')
set(H1,'LineWidth',3)
set(H2,'LineWidth',3)
axis(AX(1), [0 40 15 20]);
set(gca, 'YTick', 15:20)
axis(AX(2), [0 40 0 1]);


%
% beta function and plasma density
%
s = s_ramp;
n = n_ramp;
s = s_ft;
n = n_ft;
[AX,H1,H2] = plotyy(abscissa*100, betax_slice*100, s/1e4, n, 'plot');
set(gca, 'YTick', 1:10)
xlabel('s [cm]');
grid on;
set(H1,'Color','b')
set(H2,'Color','g')
set(get(AX(1),'Ylabel'),'String','\beta_x [cm]'); 
set(get(AX(2),'Ylabel'),'String', 'n_p / n_0 [-]');   
set(get(AX(1),'Ylabel'),'FontSize',22); 
set(get(AX(2),'Ylabel'),'FontSize',22); 
set(get(AX(1),'Ylabel'),'Color','blue'); 
set(get(AX(2),'Ylabel'),'Color','green'); 
%set(H1,'LineStyle','-')
%set(H2,'LineStyle','-')
set(H1,'LineWidth',3)
set(H2,'LineWidth',3)
axis(AX(1), [0 40 0 10]);
axis(AX(2), [0 40 0 1]);


%
% two-beam efficiency
%
do_extra_plot = 0;
if(do_extra_plot)
 N_db = qp(1).PP(1).N0;
 N_wb = qp(1).PP(2).N0;
 
  mean_E = my_flatten_struct(qp, 'mean_E');
  std_E = my_flatten_struct(qp, 'std_E');
  abscissa = [qp(:).s_timestep];

  plot(abscissa*100, mean_E(:,1), '-xr');
  xlabel('s [cm]');
  ylabel('<E> [GeV]');
  grid on;
  hold on;
  plot(abscissa*100, mean_E(:,2), '-ob');
legend('DB', 'WB', 'Location', 'Best');
%  plot(abscissa*100, mean_E(:,1)*N_db + mean_E(:,2)*N_wb, '-k'); % extract info here
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
 E_final_db = mean_E(end,1)*N_db
 E_final_wb = mean_E(end,2)*N_wb
 E_initial_db = mean_E(1,1)*N_db
 E_initial_wb = mean_E(1,2)*N_wb
 dE_db = E_initial_db-E_final_db
 dE_wb = E_final_wb-E_initial_wb
 final_trans_eff = dE_wb / dE_db * 100
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
% matched/no matched comparison
%
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




stop


%
% QP to ELE
%
%load /Users/eadli/Dropbox/SLAC/quickpic/myData/10mm_Dx0_LB_matched.mat
%pp_qp = qp_BEAMS;
load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3e17_30cm.mat
%load /Users/eadli/Dropbox/SLAC/quickpic/myData/20x20x20_match.mat
charge = 3e-9; % NB: charge is not written/transferred!
pp_qp = qp(1).PP(1).BEAM;
%[Y, I] = sort(pp_qp(:,3)); % sort according to z
%pp_qp = pp_qp(I, :);
pp_qp = my_reduce_dist(pp_qp, 1); % reduce by a factor N for quicker running
%pp_qp = pp_qp(5e4:6e4,:);
pp_qp(:,6) = 20e0 * (1 + 1e-7*randn(size(pp_qp(:,6)))); % to make monoenergetic beam
%pp_qp(:,6) = pp_qp(:,6) * 100;
filename = '/Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc';
my_write_qp2elefile(pp_qp, filename);

% analyse
%  test: re-analyse same elegant beam with quickpic
% 1) org qp beam
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam(pp_qp, [1 1 1 0 1 0]);

%
% here: propagate through elegant lens
%
% cd .../FACETSIM/ELESIM/
% sh ./track_lens.sh 
%

% 2) read in ele beam
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.out';
pp_qp = my_read_elefile2qp(filename);
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam(pp_qp, [1 1 1 0 1]);


% 4) twiss
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.twi';
[twiss, param] = my_read_elefile2twiss(filename);
%n_lensend = (64*3+2);
n_lensend = length(twiss);
ax = sqrt(twiss(end,2) / twiss(n_lensend,2))
ay = sqrt(twiss(end,4) / twiss(n_lensend,4))
%plot(twiss(1:end,1), twiss(1:end,2));
frac_plot = n_lensend / length(twiss);
%frac_plot = 6/6;
plot(twiss(1:end*frac_plot,1), twiss(1:end*frac_plot,2), 'b')
hold on;
plot(twiss(1:end*frac_plot,1), twiss(1:end*frac_plot,4), 'r')
hold off;
xlabel('s [m]');
ylabel('beta [m]');
legend('x', 'y');
grid on;

filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.mat';
% check imaging condition
[M, M_arr] = my_read_elefile2mat(filename);
m12 = M(1,2) 
m34 = M(3,4) 



%
% osiris 
%
%data = load('~/temp/test.dat');
data = load('~/temp/test2.dat'); % Feb 2012
data = load('~/temp/beforeramp.dat'); % Oct 2012
x = data(1:7:end);
y = data(2:7:end);
z = data(3:7:end);
xp = data(4:7:end);
yp = data(5:7:end);
E = data(6:7:end);
q = data(7:7:end);
% fixing offsets
%x = x - mean(x);
%y = y - mean(y);
%xp = xp - mean(x);
%yp = yp - mean(y);
pp_full = [x y z xp yp E];
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam([pp_full q], [1 1 1 0 1 0]);

% downsampling
dist_red = my_reduce_dist([pp_full q], 3, 1504);
%dist_red = [pp_full q];
pp_os = dist_red(:,1:6);
q_red = dist_red(:,7);
filename = '/Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc';
my_write_qp2elefile(pp_os, filename);
q = q_red;

% take out low charge tails and create weighted dist
if(1)
[q_sort, I] = sort(q); % sort according to q
cutfrac = 5/6;
charge_pers = sum(q_sort(1:round(end*cutfrac)) / sum(q)) * 100
n_new_particles = sum(sort(q_sort(end*cutfrac:end)) / min(q_sort(end*cutfrac:end)));
q_red_scaled = (sort(q_sort(round(end*cutfrac:end))) / min(q_sort(round(end*cutfrac:end))));
%semilogy(q_red_scaled)
pp_red = pp_os(I, :);
pp_red = pp_red(round(end*cutfrac):end, :);
%
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam(pp_red, [1 1 1 0 1 0]);
% multiply number of particle for real scaling)
pp_mult = my_mult_part(pp_red, round(q_red_scaled));
end% if

% QS1DUMP
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/ELESIM/qs1.out';
pp_qp = my_read_elefile2qp(filename);
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam(pp_qp, [1 1 1 0 1]);


% 2) read in ele beam
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.out';
pp_wafer = my_read_elefile2qp(filename);
%  sdds2stream ../ELESIM/aerogel.out  -columns=particleID  > /tmp/dat.txt
wafer = load('/tmp/dat.txt');
wafer_id = wafer(:,1);
q_wafer = q_red(wafer_id);
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam([pp_wafer q_wafer], [1 1 1 0 1]);



% lost particle id
%  sdds2stream ../FACET_WORK/facet.lost -columns=particleID,s,p  > /tmp/dat.txt
lost = load('/tmp/dat.txt');
lost_id = lost(:,1);
charge_lost = sum(q_red(lost_id)) / sum(q_red) * 100




% plot charge
semilogy(sort(q_red)/ max(q_red));
xlabel('osiris particle # (sorted for charge');
ylabel('relative charge');







%
% Li/Rb comparison
%
if(0)
%
set(gcf, 'Color', 'w');
set(0,'defaultaxesfontsize',20);
% plot some quantitites along plasma
%
%[AX, H1, H2] = plotyy([qp.s_timestep]*100, [qp.FEZ_max], [qp.s_timestep]*100, [qp.QEB_max], 'plot');
%set(H1,'Color', 'r');
%set(H2,'Color', '-b');
%set(H1,'Color', 'r');
%set(H2,'Color', '-b');
%load -mat ~/temp/l2040_ev.mat % 2040
%load -mat ~/temp/l4040_ev.mat
load -mat ~/temp/l4962_ev.mat
s_li = [qp.s_timestep];
FEZ_li = [qp.FEZ_max];
%load -mat ~/temp/r2040_ev.mat % 2040
%load -mat ~/temp/r4040_ev.mat
load -mat ~/temp/r4962_ev.mat
s_rb = [qp.s_timestep];
FEZ_rb = [qp.FEZ_max];
%load -mat ~/temp/4025match.mat
hh = plot([s_li]*100, [FEZ_li] + 0 /  (max(FEZ_li)+eps), '-xk');
set(hh, 'LineWidth', 3);
hold on;
%plot([qp.s_timestep]*100, [qp.QEB_max]/  max([qp.QEB_max]), '-xb');
%plot([qp.s_timestep]*100, [qp.QEB_max]/10, '-ok');
%load -mat ~/temp/nomatch.mat
%load -mat ~/temp/r2040_ev.mat
%qp = qp(1:end-1);
%max_field = max([qp.FEZ_max]);
hh = plot([s_rb]*100, [FEZ_rb] + 0 /  (max([FEZ_rb])+eps), ':xr');
set(hh, 'LineWidth', 3);
%plot([qp.s_timestep]*100, [qp.QEB_max]/  max([qp.QEB_max]), '-ok');
hold off;
grid on;
xlabel('s [cm]');
ylabel('Peak dec. E_z [GV/m]');
%ylabel('E_z [GV/m]  and  n_b [10 n_{p0}]');
%ylabel('E_z / E_{z,max}    ,    n_b / n_{b,max}');
%legend('Li, 20x20x40 um^3 (unmatched)', 'Rb, 20x20x40 um^3 (unmatched)');% 2040
%legend('Li, 40x40x40 um^3 (unmatched)', 'Rb, 40x40x40 um^3 (unmatched)');
legend('Li, 49x62x40 um^3 (unmatched)', 'Rb, 49x62x40 um^3 (unmatched)');
%legend('E_z, a) 25x25x40 um^3, matched', 'n_{b}, a)  25x25x40 um^3, matched');
title('PIC sims of E_z(s), for Q = 3 nC, n_0 = 2.5x10^{17} /cm^3, Field-ionized Li and Rb');
myaxis = axis;
axis([0 40 myaxis(3) myaxis(4)]);

end% if




%
% witness bunch metrics - THESE ARE THE GOOD CALCS (Nov 12, 2012)
%
datadir(26:end-1)
E0_DB = qp(1).PP(1).gamma0*0.511e-3; 
E0_WB = qp(1).PP(2).gamma0*0.511e-3; 
E1_DB = qp(end).PP(1).mean_E
E1_WB = qp(end).PP(2).mean_E
eps0 = [100.0e-6  100.0e-6];
N_dump = 2^18 / 2;
frac_part =  size(qp(1).PP(2).BEAM, 1) / N_dump 
n_qd_qw = qp(1).PP(1).N0 / qp(1).PP(2).N0 % ratio charge drive-bunch witness-bunch
E_E0 = qp(end).PP(2).mean_E / E0_WB   
sigE_E = qp(end).PP(2).std_E/ qp(end).PP(2).mean_E * 100
eps_eps0 = qp(end).PP(2).emitt ./ eps0
% energy gained versus energy extracted from DB
eff_d2w = (E1_WB - E0_WB) / (E0_DB - E1_DB) / n_qd_qw * frac_part * 100
% energy gained versus all DB energy
eff_all2w = (E1_WB - E0_WB) / (E0_DB) / n_qd_qw * frac_part * 100



%
% witness bunch analysis
%
my_ana_beam_slice(qp(end).PP(2).slice);
b = qp(1).PP(2).BEAM;
my_ana_beam(b, [1 1 1 0 1]);



%
% September 2012
% 2-bunch performance
%
met_E_E0 = [2.00  1.87   1.49  1.37   1.95  1.63   1.45   1.04  1.87 1.41 1.32 1.0 ]; % last point : 1.18, but with only 20% of particles
met_sigE_E = [6.19  5.74   5.04   4.64  5.46   11.4 5.10  7.21   3.66 7.14 4.25 5.24 ];
met_eps_eps0x = [0.99  14.1  0.99   1.64      0.99   44.4  0.99  5.61   0.99 26.4 0.99 3.04 ];
met_eps_eps0y = [0.98  19.8  0.98  1.08   0.99    49.1  0.99   5.82   0.98 37.9 0.98 6.76 ];
met_eff_d2w = [0.70 0.54  0.59  0.59   0.58  0.49   0.49    0.13  0.54   0.39    0.41  0.0]; % last point : 0.15, but with only 20% of particles

%met_eff_d2w = [0.70   0.71  0.59   0.59    0.58   0.62 0.62  0.13   0.55 0.55 0.55 0.55 ];
    
set(gcf, 'Color', 'w');
ordinate = met_E_E0;
ordinate = met_sigE_E;
ordinate = met_eff_d2w*100;
hc = plot(1:4, ordinate(1:4), '-rx');
set(hc, 'MarkerSize', 25)
hold on;
hc = plot(1:4, ordinate(5:8), '-bo');
set(hc, 'MarkerSize', 25)
hc = plot(1:4, ordinate(9:12), '-k+');
set(hc, 'MarkerSize', 25)
hold off;
legend('Pre-ionized plasma','Field-ionized Rb', 'Field-ionized Li');
%xlabel('scenario');
title('Energy gain after 70 cm of plasma at n_0 = 1x10^{17}/cm^3, E_0 = 25 GeV'); ylabel('Witness bunch <E> / <E_0> [-]');  myaxis = axis; axis([myaxis(1) myaxis(2) 1 2]);
title('Energy spread after 70 cm of plasma at n_0 = 1x10^{17}/cm^3, E_0 = 25 GeV'); ylabel('Witness bunch \sigma_E / E [%]');  myaxis = axis; axis([myaxis(1) myaxis(2) 0 15]);
title('Energy efficiency after 70 cm of plasma at n_0 = 1x10^{17}/cm^3, E_0 = 25 GeV'); ylabel('Drive to witness bunch efficieny [%]'); myaxis = axis; axis([myaxis(1) myaxis(2) 0 100]);
%set(0,'DefaultTextInterpreter', 'latex')
set(gca,'XTick',[1:1:4]);
set(gca, 'xticklabel', {'Q=6nC,r=3.3um(mat)'; 'Q=6nC,r=20um(unmat)'; 'Q=3nC,r=3.3um(mat)'; 'Q=3nC,r=20um(unmat)  '});





%
%
%

% Esarey example E-157
gamma = 30/.511e-3;
kb = 2.6e3 / sqrt(2*gamma); % [1e17/cm3]
r = 100e-6;
% PWFA-LC
gamma = 25/.511e-3;
gamma = 500/.511e-3;
kb = 6e4 / sqrt(2*gamma); % [1e17/cm3]
r = 0.5e-3;
r = 10e-6;
betatron_rad_loss = 0.48*(gamma/1e4)^4 * (kb/1e2)^4 * (r/1e-3)^2 % [GeV/m]

% other formula
gamma = 25/.511e-3;
kp = 6e4; % [1e17/cm3]
r = 10e-6;
betatron_rad_loss = 12.00*(gamma/1e4)^2 * (kp/1e4)^4 * (r/1e-2)^2 % [GeV/m]









%  Special case of E-profile for Rb vs Li plots.  (was in main loop
%  of QuickPIC ana, but taken away since it was not general) :
%
%  plot Eprof
%
if(do_Eprof)
  % assume Rb pre-saved
% save -mat ~/temp/final_Rb.mat n E n_dec_frac n_acc_frac
  load('~/temp/final_Rb');
  n_Rb = n;
  E_Rb = E;
  n_dec_frac_Rb = n_dec_frac;
  n_acc_frac_Rb = n_acc_frac;

  % assume Li current
  n_hist = 101;
  [n,E] = hist(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6), n_hist);
  % calc acc and dec charge
  sigma_E0 = 0.04; % 1 sigma 
  % to compare with experiment: take 2 sigma away from center
  E0 = qp(n_3D_counter).gamma0*0.511e-03;
  E_dec_max = E0*(1-2*sigma_E0);
  E_acc_min = E0*(1+2*sigma_E0);
  n_dec_frac = sum( n(1:min(find(E>E_dec_max))) ) / sum(n);
  n_acc_frac = sum( n(min(find(E>E_acc_min)):end) ) / sum(n);
  
  subplot(2,3,3);
  hh = semilogy([E0 E0+eps], [1 1e4], '-b');  
  set(hh,'LineWidth',5)
  hold on;
  hh = plot([E_dec_max E_dec_max+eps], [1 1e4], '-g');  
  set(hh,'LineWidth',2)
  hh = plot([E_acc_min E_acc_min+eps], [1 1e4], '-g');  
  set(hh,'LineWidth',2)
  % Li
  hh = semilogy(E,n / max(max(n), max(n_Rb)) * 1e4, '-k');
  set(hh,'LineWidth',3)
  % Rb
  hh = semilogy(E_Rb,n_Rb / max(max(n), max(n_Rb)) * 1e4, '-r');
  set(hh,'LineWidth',3)
  hold off;
  grid on;
  xlabel('E [GeV]');
  ylabel('count [arb.units.]');
  axis([10 40 1 1e4]);
  % 
  set(0,'defaulttextfontsize',18);
%  set(0,'defaulttextfontname','Courier');
  text(31, 100, ['Li-%: ' sprintf('%0.2g', n_dec_frac*100)]);
  text(31, 30, ['Li+%: ' sprintf('%0.2g', n_acc_frac*100)]);
  text(31, 10, ['Rb-%: ' sprintf('%0.2g', n_dec_frac_Rb*100)]);
  text(31, 3, ['Rb+%: ' sprintf('%0.2g', n_acc_frac_Rb*100)]);
  hh = legend('Init. beam', 'Dec/Acc. cut', '', 'Final Li', 'Final Rb');
  set(hh,'FontSize',12);
%  set(hh,'FontName','Courier');
  
end% if plot


%
%
%
load ~/temp/beforeramp.dat;
B = reshape(beforeramp, length(beforeramp)/7, 7);
 
plot(B(:,3)/1e3, B(:,6), 'x');
xlabel('z [mm]');
xlabel('E [GeV]');




%
%
%
emitt_x = my_flatten_struct(qp, 'emitt(1)');
emitt_y = my_flatten_struct(qp, 'emitt(2)');
abscissa = [qp(:).s_timestep];
plot(abscissa*1e2, emitt_x*1e6, '-r');
hold on;
plot(abscissa*1e2, emitt_y*1e6, '-b');
hold off;
grid on;
xlabel('s [cm]');
ylabel('Total \epsilon_N [um]');
legend('x', 'y');




%
% Summary charge scan, Nov 13, 2012
%
met_E_E0 = [1.07 1.17 1.40 1.68];
met_sigE_E = [1.8 4.4 8.8 16.7];
met_eff_d2w = [36 43  45 47.8];
%
PI_E_E0 = [1.50 1.59 1.59 1.69];
PI_sigE_E = [9.8 12 12 20];
PI_eff_d2w = [45 45 48 48];

%
% FI 
%
set(gcf, 'Color', 'w');
ordinate = met_E_E0;
%ordinate = met_sigE_E;
%ordinate = met_eff_d2w;
hc = plot(1:4, ordinate(1:4), '-rx');
set(hc, 'MarkerSize', 25)
set(gca,'XTick',[1:1:4]);
set(gca, 'xticklabel', {'Q=1.5x10^10e'; 'Q=2.0x10^10e'; 'Q=3.0x10^10e'; 'Q=4.0x10^10e  '});
ylabel('Witness bunch <E> / E_0');
title('Energy gain vs. total charge (after coll.), Field-ionized Rb, E_0 = 20 GeV, \sigma_r = 20 um'); % ylabel('Witness bunch <E> / <E_0> [-]');  myaxis = axis; axis([myaxis(1) myaxis(2) 1 2]);
grid on;

hold on;
%
% FI + PI
%
PI_E_E0 = [1.50 1.59 1.59 1.69];
PI_sigE_E = [9.8 12 12 20];
PI_eff_d2w = [45 45 48 48];

set(gcf, 'Color', 'w');
hold on;
ordinate = met_E_E0;
hc = plot(2:5, ordinate(1:4), '-rx');
set(hc, 'MarkerSize', 25)
hold on;
ordinate = PI_E_E0;
hc = plot(1:1, ordinate([1 ]), '-ob');
set(hc, 'MarkerSize', 25)
hc = plot(1:2, ordinate([2 4]), '-xm');
set(hc, 'MarkerSize', 25)
hold off
set(hc, 'MarkerSize', 25)
set(gca,'XTick',[1:1:5]);
set(gca, 'xticklabel', {'Q=1.2x10^10e'; 'Q=1.5x10^10e'; 'Q=2.0x10^10e'; 'Q=3.0x10^10e'; 'Q=4.0x10^10e  '});
axis([1 5 1 1.8]);
ylabel('Witness bunch <E> / E_0');
title('Energy gain vs. total charge (after coll.), E_0 = 20 GeV, \sigma_r = 20 um'); % ylabel('Witness bunch <E> / <E_0> [-]');  myaxis = axis; axis([myaxis(1) myaxis(2) 1 2]);
grid on;
legend('Field-ion Rb, 300/30', 'Pre-ion, 300/30', 'Pre-ion, 100/100', 'Location', 'Best');



%
% calc dist
%
pp = qp(1).PP(2).BEAM;
[N,X ] = hist(pp(:,6), 201); 
semilogy(X, N/max(N));
xlabel('E [GeV]'); 
ylabel('count [a.u.]');
grid on;

E_mean = mean(pp(:,6));
frac = 0.03;
n_above = abs(X - E_mean) < (E_mean*frac)
N_above = N(n_above);
sum(N_above) / sum(N)
xlabel('E [GeV]');
 
hold on;
 hh = plot([mean_E+E_mean*frac mean_E+E_mean*frac+eps], [1e-5 1e0], '-g');  
 set(hh,'LineWidth',2)
 hh = plot([mean_E-E_mean*frac mean_E-E_mean*frac+eps], [1e-5 1e0], '-g');  
 set(hh,'LineWidth',2)
hold off;



%
%
%
plot([qp(:).s_timestep]*100, my_flatten_struct(qp, 'x_head'), '-ob');
hold on;
plot([qp(:).s_timestep]*100, my_flatten_struct(qp, 'x_ion_head'), '-xr');
plot([qp(:).s_timestep]*100, my_flatten_struct(qp, 'z_ion_head'), '-+g');
hold off;
grid on;
xlabel('s [cm]');
ylabel('<x>, z [um]');
legend('Eroded beam','New head (not eroded)',  'z new head', 'Location', 'Best');
title_text = datadir;
title(title_text, 'interpreter', 'None');











% at 40 cm
z0_0 = -52.0; % assumed start of beam w/o head erosion
%z0_0 = -45.2; % assumed start of beam w/o head erosion
% FI emittance matched
FI_z040 = z0_0 - fliplr([5.2 -15.6 -30.4 -38.2 -41.6 -44.2 -44.2 -44.2 NaN NaN]);
FI_ctm40 = fliplr([-96.4 -63.9 -42.6 -30.22 -22.2 -16.4 -8.9 -4.7 NaN NaN]);
% PI emittance matched
FI_z040 = z0_0 - fliplr([-27.8 -37.3 -42.5 -46.0 -46.0  -47.7 -51.2 NaN   NaN  -52.0]);
FI_ctm40 = fliplr([-26.6 -21.4 -16.4 -10.7  -8.4  -6.7  -5.6  NaN NaN -5.6]);
% FI, constant sr=5
FI_z040 = z0_0 - fliplr([ 5.2  -13.9  -19.1 -35.6  -33.8 -33.8  NaN -32.9  NaN NaN]);
FI_ctm40 = fliplr([ -96.4 -32.8 -14.2 -7.0 -4.6 -3.8  NaN -3.1 NaN NaN]);
% PI, constant sr=5
%FI_z040 = z0_0 - fliplr([-27.8 -30.4 -41.6  -39.0 -39.0  NaN  NaN -39.0 NaN -39.0 ]);
%FI_ctm40 =  fliplr([-26.6 -13.3  -7.1 -4.8 -4.3  NaN  NaN -4.2 NaN  -4.2]);

sr = 5./sqrt([2^99 2^18 2^14 2^12 2^8 2^5 2^4 2^3 2^2 2^1 2^0]);
sr = sr(2:end);
%FI_z040 = [-1 FI_z040]; 
ab_pow = 4/2;
or_pow = 2/2;
abscissa = sr.^(ab_pow); 
n_start = 2; % 1: include art. zero, 2: skip art. zero
hh = plot(abscissa(n_start:end), (-FI_ctm40(n_start:end)).^(or_pow), '-xr');
set(hh, 'MarkerSize', 25)
hold on;
hh = plot(abscissa(n_start:end), (-FI_z040(n_start:end)).^(or_pow), '-ob');
set(hh, 'MarkerSize', 25)
hh = plot(abscissa(n_start:end), FI_ctm40(n_start:end) ./ FI_z040(n_start:end)*1, '-+m');
set(hh, 'MarkerSize', 25)
hold off;
grid on;
xlabel(['\sigma_r [um]' ' pow ' sprintf('%0.1f', ab_pow)]);
%xlabel(['prop.to. \epsilon [um]' ' pow ' sprintf('%0.1f', ab_pow.^2)]);
ylabel('Ers. front, and ctm at ers. front @ 40 cm [um]');
title('FI, emittance matched');
title('FI, emittance, \sigma_r = 5 um');
%title('PI, emittance matched');
%title('PI, emittance, \sigma_r = 5 um');
legend('ctm [um]', 'ers.front [um]', 'ctm/ers.front', 'Location', 'Best');
myaxis = axis;
%axis([0 myaxis(2) myaxis(3) 100]);
