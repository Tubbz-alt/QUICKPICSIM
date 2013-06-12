clear all;

%
% SAVE QS VALUES
%
%QS_setting = -20.35;
QS_setting = +8.0;
E0 = 20.35;

%
% construct test beams
%
% x, xp
varix = 1; % 1:x, 2:y
E0 = E0; %+1e-5;
N = 1e4;
N2 = sqrt(N);
x_max = 3e3; % [um]
xx = linspace(-x_max,x_max, N2);
varname = {'x', 'y'};
varu = [1 2];
varup = [4 5];
varun = [2 1];  
varupn = [5 4];

% use z to plot z-E
z_max = 1e-5; % [um]
zz = linspace(-z_max,z_max, N2);
dE_max = 5e-6; % [-]
dEE = linspace(-dE_max,dE_max, N2);
y_max = 5e-3; % [um]
yy = linspace(-y_max,y_max, N2);

for nx=0:(N2-1),
  for kx=0:(N2-1),
%    pp_test(kx+nx*N2+1, varu(varix)) = xx(kx+1);
    pp_test(kx+nx*N2+1, 3) = zz(kx+1);
%    pp_test(kx*N2+nx+1, varup(varix)) = xx(kx+1);
    pp_test(kx*N2+nx+1, 6) = E0*(1+dEE(kx+1));
%    pp_test(kx*N2+nx+1, 6) = E0*(1+0.01);
%    pp_test(kx*N2+nx+1, 5) = yy(kx+1);
%    pp_test(kx*N2+nx+1, 2) = yy(kx+1);
% for pure cut
    pp_test(kx+nx*N2+1, 1) = xx(kx+1)/1e6;
    pp_test(kx*N2+nx+1, 2) = xx(kx+1)/1e6;
    pp_test(kx+nx*N2+1, 4) = xx(kx+1);
    pp_test(kx*N2+nx+1, 5) = xx(kx+1);
  end% for
end% for
%pp_test(:, 6) = E0 +randn(size(pp_test,1),1)*1e-6*E0;
%pp_test(:, 3) = linspace(-z_max, z_max, size(pp_test, 1));
%pp_test(:, varun(varix)) = pp_test(:, varu(varix))/1e4;
%pp_test(:, varupn(varix)) = pp_test(:, varup(varix))/1e4;

% for E cal
%E_range = E0*(0.1:0.001:3);
%pp_test = [];
%pp_test(:, 6) = E_range;

% for cut
%u_range = 5e4:0.1:5e4
%pp_test = [];
%pp_test(:, 6) = E_range;


% temp for test
%pp_test = [];
%pp_test(:,5) = yy;
%pp_test(:,6) = E0;

%load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3e17_30cm.mat
%load /Users/eadli/Dropbox/temp/qp_hose.mat;  % hose @ 40 cm
%load /Users/eadli/Dropbox/temp/qp_hose0.mat; % initial beam
%load /Users/eadli/Dropbox/temp/qp_hose_new.mat;  %
%load /Users/eadli/Dropbox/temp/qp_hose_new_se_t0.mat;  %
load /Users/eadli/Dropbox/temp/qp_hose_new_se_t0_long.mat;  %
charge = 3e-9; % NB: charge is not written/transferred!
pp_qp = qp(1).PP(1).BEAM;
pp_qp = my_reduce_dist(pp_qp, 4, 150475); % reduce by a factor N for quicker running
%pp_qp(:,6) = pp_qp(:,6)+0; % eq to -QSBEND
% reduce emittance artifactly by factor
% RAMP, factor N in beta transfer
beta_ramp = 100;
pp_acc_indx = pp_qp(:,6) > E0;
pp_qp(pp_acc_indx,1) = pp_qp(pp_acc_indx,1) * beta_ramp;
pp_qp(pp_acc_indx,2) = pp_qp(pp_acc_indx,2) * beta_ramp;
pp_qp(pp_acc_indx,4) = pp_qp(pp_acc_indx,4) / beta_ramp;
pp_qp(pp_acc_indx,5) = pp_qp(pp_acc_indx,5) / beta_ramp;
% multiply particle count in acc tail by 3
%pp_acc = pp_qp(pp_qp(:,6) > E0, :);
%pp_qp = [pp_qp; pp_acc]; %; pp_acc; pp_acc];
%  USE LOADED BEAM BEAM INSTEAD
pp_test = pp_qp;
% ADD BEAMLET
x_div = -1.0e3; % [urad]
y_div = 1.0e3; % [urad]
dE = 0;
pp_beamlet = [0 0 0 x_div y_div E0+dE];
N_beamlet = round(size(pp_test,1)/100);
for kx=1:N_beamlet,
  beamlet_emit = 1/50;
  pp_new = [x_div/10*randn y_div*beamlet_emit*randn 0  x_div+x_div*beamlet_emit*randn y_div+y_div*beamlet_emit*randn E0+dE+randn*(E0+dE)*beamlet_emit];
  pp_beamlet = [pp_beamlet; pp_new];
end% for
pp_test = [pp_test; pp_beamlet];
%pp_test = [pp_beamlet];

filename = '/Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc';
my_write_qp2elefile(pp_test, filename);
%[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp_test, [1 1 1 0 1 0 0]);

% write QS param file
KQS1_0 = 3.077225846087095e-01;
KQS2_0 = -2.337527121004531e-01;
KQS1 = KQS1_0*(E0+QS_setting)/E0;
KQS2 = KQS2_0*(E0+QS_setting)/E0;
filename = '/tmp/QS_elegant.par.txt';
fid = fopen(filename, 'w');
fprintf(fid, 'QS1 K1 %.12e\n', KQS1);
fprintf(fid, 'QS2 K1 %.12e\n', KQS2);
fclose(fid);
system('cd /Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_OPTICS/; /Users/eadli/Dropbox/SLAC/elegant/bin/plaindata2sdds /tmp/QS_elegant.par.txt /Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_OPTICS/QS_elegant.par -col=ElementName,string -col=ElementParameter,string  -col=ParameterValue,double -noRowCount');

disp(['QS_setting: ' num2str(QS_setting) '']);
disp('Run elegant now');
% ./sim_elegant.sh -i /Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc -d ../FACET_OPTICS/facet_v27.4.QS.ele   -b FACET -n off   -r ../FACET_OPTICS/R56Params/10.0mmR56.par   -l ../FACET_OPTICS/facet_simple_dump_only.lte    -o ../FACET_WORK/facet_lens.out

pause;
%
% 4) track beam in dump line
%
% system(['cd /Users/eadli/Dropbox/SLAC/elegant/FACETSIM/ELESIM;  ./sim_elegant.sh -i /Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc -d ../FACET_OPTICS/facet_v27.4.dynamic.ele   -b FACET -n off   -r ../FACET_OPTICS/R56Params/10.0mmR56.par   -l ../FACET_OPTICS/facet_simple_dump_only.lte    -o ../FACET_WORK/facet.out']);



% analyze final beam
%
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.out';
[pp_DUMP, pp_DUMP7] = my_read_elefile2qp(filename);
if(size(pp_DUMP,1) > 0)
  [N, I] = sort(pp_DUMP7(:,7));
end% if
%pp_DUMP = pp_DUMP(I,:); % re-sort pp_DUMP on PID
%[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp_DUMP, [0 0 0  0 1 0 0]); sigx, sigy
%[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp_DUMP, [1 1 1  0 1 0 0]); sigx, sigy

my_ana_beam(pp_DUMP, [1 0 0 0  1 0 0 1]);
pause;
subplot(2,3,5);
axis([-40 40 -100 100]);
subplot(2,3,6);
axis([-40 40 0 50]);
subplot(2,3,2);
axis([-40 40 0 50]);
subplot(2,3,3);
axis([-0 50 -50 50]);
pause;

subplot(1,1,1);
x_red = 1e3; 
%plot(pp_test(:,varu(varix))/x_red, pp_test(:,varup(varix))/x_red, 'rx');
plot(pp_test(:,4)/x_red, pp_test(:,5)/x_red, 'rx'); % xp, yp
plot(pp_test(:,6), pp_test(:,5)/x_red, 'rx'); % E, yp
hold on;
if(size(pp_DUMP,1) > 0)
%  plot(pp_test(pp_DUMP7(:,7),varu(varix))/x_red, pp_test(pp_DUMP7(:,7),varup(varix))/x_red, 'kx');
%  plot(pp_test(pp_DUMP7(:,7),4)/x_red, pp_test(pp_DUMP7(:,7),5)/x_red, 'ko');
  plot(pp_test(pp_DUMP7(:,7),6), pp_test(pp_DUMP7(:,7),5)/x_red, 'ko');
end% if
%plot([-x_max x_max]/x_red, [0 eps], '-k');
%plot([0 eps], [-x_max x_max]/x_red, '-k');
hold off;
pbaspect([1 1 1]);
%xlabel(['' char(varname(varix)) ' [mm]']);
%ylabel(['' char(varname(varix)) 'p [mrad]']);
%xlabel([' xp [mrad]']);
%ylabel([' yp [mrad]']);
xlabel([' E [GeV]']);
ylabel([' yp [mrad]']);
%title(['               ' 'Aperture ' char(varname(varix)) '.  E_{part} = ' num2str(E0, 4) ' GeV, Imaging E = ' num2str(QS_setting+E0, 4) ' GeV   ("QS = ' num2str(QS_setting,4) '")']);
title(['               ' 'Aperture.  E_{part} = ' num2str(E0, 4) ' GeV, Imaging E = ' num2str(QS_setting+E0, 4) ' GeV   ("QS = ' num2str(QS_setting,4) '")']);
title(['               ' 'Aperture.  Imaging E = ' num2str(QS_setting+E0, 4) ' GeV   ("QS = ' num2str(QS_setting,4) '")']);
%title(['               ' 'Aperture. Freely propagating beam (QS quadrupoles off).']);

pause;

% TEM
my_ana_beam(pp_DUMP, [1 0 0 0  1 0 0 1]);
% get energy scale
load ~/Dropbox/temp/yE_nom.mat;

pause;

% simulated cher image
x_var = 1;
y_var = 2;
xx= (-100:0.5:100)*1e3; 
yy = (-100:0.5:40)*1e3;
dist = hist2(pp_DUMP(:,x_var), pp_DUMP(:,y_var), xx, yy);
EE = my_2d_mapping(E_nom, y_nom, yy);
%pcolor(xx/1e3, yy/1e3, log10(dist)/1e3); 
pcolor(xx/1e3, EE, log10(dist)/1e3); 
shading('flat');  
colorbar;
mycaxis=caxis; 
caxis([mycaxis(1) mycaxis(2)/1]);
axis([-50 50 10 30])
xlabel('x_{sim} [mm]');
ylabel('E_{meas,sim} [GeV]');

pause;

% get energy scale
load ~/Dropbox/temp/yE_nom.mat;
pp_DUMP7(:,8) = my_2d_mapping(E_nom, y_nom, pp_DUMP7(:,2)); % map sim y-axis to sim E-axis 
%  plot(pp_test(pp_DUMP7(:,7),2)/x_red, pp_test(pp_DUMP7(:,7),6)/x_red, 'bx');
n_xscan = 1:N2;
n_xpscan = 1:N2:length(pp_test);
%plot(pp_DUMP7(n_xscan,2)/x_red, pp_DUMP7(n_xscan,6), 'kx');
plot(pp_DUMP7(n_xpscan,2)/x_red, pp_DUMP7(n_xpscan,6), 'kx');
xlabel('y_{CHER} [mm]');
ylabel('E [GeV]');
%plot(pp_test(n_xpscan,2)/x_red, pp_DUMP7(n_xpscan,2)/x_red, 'kx');
%ylabel('y_{CHER} [mm]');
%xlabel('y_{OP} [mm]');
%plot(pp_test(n_xpscan,5)/x_red, pp_DUMP7(n_xpscan,2)/x_red, 'kx');
%plot(pp_test(:,5)/x_red, pp_DUMP7(:,2)/x_red, 'kx');
%ylabel('y_{CHER} [mm]');
%xlabel('yp_{OP} [mrad]');
%plot(pp_test(n_xpscan,6)/x_red, pp_DUMP7(n_xpscan,2)/x_red, 'kx');
plot(pp_DUMP7(n_xpscan,8), pp_DUMP7(n_xpscan,2)/x_red, 'kx');
ylabel('E_{CHER} [GeV]');
xlabel('yp_{OP} [mrad]');
grid on

% interpolate energy
%y_nom = pp_DUMP(:,2);
%E_nom = pp_DUMP(:,6);
%save ~/Dropbox/temp/yE_nom.mat y_nom E_nom;


stop






%
% CHECK IMAGING CONDITION
%
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.mat';
[M, M_arr] = my_read_elefile2mat(filename);
m12 = M(1,2) 
m34 = M(3,4) 
ax = sqrt(twiss(end,2) / twiss(1,2))
ay = sqrt(twiss(end,4) / twiss(1,4))

%
% for analysis: twiss in plasma lens
%

filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.twi';
[twiss] = my_read_elefile2twiss(filename);
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.mat';
[M, M_arr] = my_read_elefile2mat(filename);
twiss(1,2) , twiss(1,4)
twiss(end,2); , twiss(end,4);
n_lensend = length(twiss);
m12 = M(1,2) 
m34 = M(3,4) 
ax = sqrt(twiss(end,2) / twiss(1,2))
ay = sqrt(twiss(end,4) / twiss(1,4))
%frac_plot = n_lensend / length(twiss);
frac_plot = 9/10;
frac_plot = 1;
plot(twiss(1:end*frac_plot,1), twiss(1:end*frac_plot,2), 'b')
hold on;
plot(twiss(1:end*frac_plot,1), twiss(1:end*frac_plot,4), 'r')
hold off;
xlabel('s [m]');
ylabel('beta [m]');
legend('x', 'y');
grid on;



%
% Calc ax, ay
%
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.twi';
[twiss, param] = my_read_elefile2twiss(filename);
%n_lensend = (64*3+2);
n_lensend = 1; %length(twiss);
%n_lensend = length(twiss);
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


%
%
%
%    p0 = 4 GeV
%E = [1 2
%n_loss = [38032  0 
%    p0 = 8 GeV
%E = [2 3
%n_loss = [19801  0 
%    p0 = 10 GeV
%E = [3 4
%n_loss = [186 0
%    p0 = 20 GeV
%E = [4 5 6
%n_loss = [97120 4038 0
     