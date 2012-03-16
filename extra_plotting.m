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
pp_qp = my_reduce_dist(pp_qp, 2); % reduce by a factor N for quicker running
%pp_qp = pp_qp(5e4:6e4,:);
%pp_qp(:,6) = 20e0 * (1 + 1e-6*randn(size(pp_qp(:,6)))); % to make monoenergetic beam
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
data = load('~/temp/test2.dat');
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
dist_red = my_reduce_dist([pp_full q], 1, 1504);
%dist_red = [pp_full q];
pp_os = dist_red(:,1:6);
q_red = dist_red(:,7);
filename = '/Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc';
my_write_qp2elefile(pp_os, filename);


% take out low charge tails
if(0)
[q_sort, I] = sort(q); % sort according to q
cutfrac = 2/6;
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