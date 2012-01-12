clear all;
clf;

%
% FACETSIM ana :
%
%
% update working dir to your QUICKPICSIM folder
working_dir = '/Users/eadli/Dropbox/SLAC/quickpic/QUICKPICSIM/'; eval(['run ' working_dir 'my_SI_params.m']); % import my standard SI constants

% update elegant exe dir to your elegant exe dir
elexedir = '/Users/eadli/Dropbox/SLAC/elegant/bin/';

% update elegant data output dir to your elegant data output dir
myeledir = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/';
hz_Nstd_init = 3;

rpinput_template_file = [working_dir 'rpinput/rpinput_template'];
rpinput_output_file =  [working_dir 'rpinput/rpinput'];


%
% user define whether 1 or 2 QP bunches should be simulated 
% (this will split the bunch dist at local min)
%
N_bunch_qp = 1;

%
%
%
do_pause = 1;


%
% temp
%
% "NEW" : two bunch collim params (scattered all over)
%myeledir = '/Users/eadli/Dropbox/SLAC/elegant/elegant_dumps/facet_new/'
%hz_Nstd_init = 0.3;
%N_bunch_qp = 2;

% "ORIG" : two bunch collim params (even worse)
%myeledir = '/Users/eadli/Dropbox/SLAC/elegant/elegant_dumps/facet_orig/'
%hz_Nstd_init = 0.03;
%N_bunch_qp = 2;



%
% load data
%
myfile1 = [myeledir 'facet.out'];
%myfile1 = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_BEAMS/facet.in';
myfile2 = [myeledir 'part.1.1.dat'];
system([elexedir 'sdds2stream ' myfile1 ' -col=x,xp,y,yp,t,p,particleID > ' myfile2]);
pp = load(myfile2);
myfile1 = [myeledir 'facet.out'];
myfile2 = [myeledir 'charge.1.1.dat'];
system([elexedir 'sdds2stream ' myfile1 ' -par=Charge > ' myfile2]);
charge = load(myfile2);
myfile1 = [myeledir 'facet.twi'];
myfile2 = [myeledir 'twiss.1.1.dat'];
system([elexedir 'sdds2stream ' myfile1 ' -col=s,betax,alphax,psix,etax,etaxp,xAperture,betay,alphay,psiy,etay,etayp,yAperture > ' myfile2]);
tw = load(myfile2);
% converting t to z
pp(:,5) = pp(:,5) * SI_c;
% taking out mean z
pp(:,5) =  (pp(:,5)-mean(pp(:,5)));
% converting p [m_e*c] to E [ GeV ] (assuming ultrarel)
pp(:,6) = pp(:,6) * (SI_em*SI_c^2/SI_e);
% define r = sqrt(x^2 + y^2) 
rr = sqrt(pp(:,1).^2 + pp(:,3).^2);

% plot twiss beta
plot(tw(:,1), tw(:,2), 'b');
hold on;
plot(tw(:,1), tw(:,8), 'r');
hold off;
xlabel('s [m]');
ylabel('\beta [m]');
legend('\beta_x', '\beta_y');
grid on;
% plot twiss eta
plot(tw(:,1), tw(:,5), 'b');
hold on;
plot(tw(:,1), tw(:,11), 'r');
hold off;
xlabel('s [m]');
ylabel('\eta [m]');
legend('\eta_x', '\eta_y');
grid on;


% final dispersion [mm]
Dx = tw(end,5) * 1e3
betax = tw(end,2);
alphax = tw(end,3);
Dy = tw(end,11) * 1e3
betay = tw(end,8);
alphay = tw(end,9);





% hist params
%
hz_N = 51;
hz_Nstd = hz_Nstd_init;
% one bunch params
%hz_Nstd = 3;
%hz_N = 31;

% hist limits
hz_mean = mean(pp(:,5));
hz_std = std(pp(:,5));
z_min = (hz_mean-hz_Nstd*hz_std);
z_max = (hz_mean+hz_Nstd*hz_std);
hist_z = z_min:(2*hz_Nstd*hz_std)/hz_N:z_max;
[NNz_total, XXz_total] = hist(pp(:,5), hist_z);
% disregard outliers
NNz_total(1) = 0;
NNz_total(end) = 0;
% visualize full dist
subplot(2,2,1);
my_densplot(my_ele2qp(pp), 3, 6, 0, hz_Nstd);
subplot(2,2,3);
bar(XXz_total*1e6, NNz_total);
xlabel('z [um]'); 
ylabel('count [a.u.]');
myaxis = axis;
axis([z_min*1e6  z_max*1e6  myaxis(3)  myaxis(4)]);
grid on;

if( do_pause )
  pause;
end% if


if(N_bunch_qp == 2)
  % find local minimum in order to split bunch into two QP bunches
  n_min = my_find_local_min(NNz_total)
  z_min = XXz_total(n_min)
  [pp_split(1).pp, pp_split(2).pp, N1, N2] = my_split_dist(pp, z_min);
  charge_tot = charge;
  charge(1) = charge_tot * N1 / (N1+N2);
  charge(2) = charge_tot * N2 / (N1+N2);
end% if

for n_bunch_qp=1:N_bunch_qp,
  if(N_bunch_qp == 2)
    pp = pp_split(n_bunch_qp).pp;
    
    % redo Z hist for each bunch
    if(n_bunch_qp == 2)
      % second bunch has all these outliers
      hz_Nstd = hz_Nstd_init;
      hz_N = 31;
    else
      hz_Nstd = 3;
      hz_N = 31;
    end% if
    hz_mean = mean(pp(:,5))
    hz_std = std(pp(:,5))
    z_min = (hz_mean-hz_Nstd*hz_std)
    z_max = (hz_mean+hz_Nstd*hz_std)
    hist_z = z_min:(2*hz_Nstd*hz_std)/hz_N:z_max;
    [NNz_bunch, XXz_bunch] = hist(pp(:,5), hist_z);
    NNz(n_bunch_qp, :) = NNz_bunch;
    XXz(n_bunch_qp, :) = XXz_bunch;
    % disregard outliers
    NNz(n_bunch_qp, 1) = 0;
    NNz(n_bunch_qp, end) = 0;
  else
    %   simply keep original bunch
    % pp = pp; 
    NNz = NNz_total;
    XXz = XXz_total;
  end% if
  
hx_mean = mean(pp(:,1));
hx_std = std(pp(:,1));
hx_Nstd = hz_Nstd;;
hx_Nstd = hz_Nstd;
if(hx_Nstd < 0.3 )
  hx_Nstd = 0.3;
end% if
hx_N = hz_N;
x_min = (hx_mean-hx_Nstd*hx_std);
x_max = (hx_mean+hx_Nstd*hx_std);
hist_x = x_min:(2*hx_Nstd*hx_std)/hx_N:x_max;
[NNx, XXx] = hist(pp(:,1), hist_x);
% disregard outliers
NNx(1) = 0;
NNx(end) = 0;

hy_mean = mean(pp(:,3));
hy_std = std(pp(:,3));
hy_Nstd = hz_Nstd;
hy_N = hz_N;
y_min = (hy_mean-hy_Nstd*hy_std);
y_max = (hy_mean+hy_Nstd*hy_std);
hist_y = y_min:(2*hy_Nstd*hy_std)/hy_N:y_max;
[NNy, XXy] = hist(pp(:,3), hist_y);
% disregard outliers
NNy(1) = 0;
NNy(end) = 0;


hxp_mean = mean(pp(:,2));
hxp_std = std(pp(:,2));
hxp_Nstd = hz_Nstd;
hxp_N = hz_N;
xp_min = (hxp_mean-hxp_Nstd*hxp_std);
xp_max = (hxp_mean+hxp_Nstd*hxp_std);
hist_xp = xp_min:(2*hxp_Nstd*hxp_std)/hxp_N:xp_max;
[NNxp, XXxp] = hist(pp(:,2), hist_xp);
% disregard outliers
NNxp(1) = 0;
NNxp(end) = 0;

hyp_mean = mean(pp(:,4));
hyp_std = std(pp(:,4));
hyp_Nstd = hz_Nstd;
hyp_N = hz_N;
yp_min = (hyp_mean-hyp_Nstd*hyp_std);
yp_max = (hyp_mean+hyp_Nstd*hyp_std);
hist_yp = yp_min:(2*hyp_Nstd*hyp_std)/hyp_N:yp_max;
[NNyp, XXyp] = hist(pp(:,4), hist_yp);
% disregard outliers
NNyp(1) = 0;
NNyp(end) = 0;


subplot(2,3,1);
my_densplot(my_ele2qp(pp), 3, 6, 0, hz_Nstd);

subplot(2,3,3);
% fit gaussian to get std deviation of a gaussian distribution
fit_X = XXz(n_bunch_qp, :);
fit_Y = NNz(n_bunch_qp, :);
save -mat /tmp/fitfunc.dat fit_X fit_Y
init_guess = [1e3 0 100e-6];
result = fminsearch('mygaussfit2', init_guess);
gauss_sigma_z(n_bunch_qp) = result(3)
gauss_mu_z = result(2)
gauss_A_z = result(1)/gauss_sigma_z(n_bunch_qp)/sqrt(2*pi);
xlabel('z [m]'); 

subplot(2,3,2);
fit_X = XXx;
fit_Y = NNx;
save -mat ~/Dropbox/temp/fitfunc.dat fit_X fit_Y
init_guess = [1e3 0 10e-6];
result = fminsearch('mygaussfit2', init_guess);
gauss_sigma_x(n_bunch_qp) = result(3)
gauss_mu_x = result(2)
gauss_A_x = result(1)/gauss_sigma_x(n_bunch_qp)/sqrt(2*pi);
xlabel('x [m]'); 

subplot(2,3,4);
fit_X = XXy;
fit_Y = NNy;
save -mat ~/Dropbox/temp/fitfunc.dat fit_X fit_Y
init_guess = [1e3 0 10e-6];
result = fminsearch('mygaussfit2', init_guess);
gauss_sigma_y(n_bunch_qp) = result(3)
gauss_mu_y = result(2)
gauss_A_y = result(1)/gauss_sigma_y(n_bunch_qp)/sqrt(2*pi);
xlabel('y [m]'); 


subplot(2,3,5);
fit_X = XXxp;
fit_Y = NNxp;
save -mat ~/Dropbox/temp/fitfunc.dat fit_X fit_Y
init_guess = [1e3 0 10e-5];
result = fminsearch('mygaussfit2', init_guess);
gauss_sigma_xp(n_bunch_qp) = result(3)
gauss_mu_xp = result(2)
gauss_A_xp = result(1)/gauss_sigma_xp(n_bunch_qp)/sqrt(2*pi);
xlabel('xp [rad]'); 


subplot(2,3,6);
fit_X = XXyp;
fit_Y = NNyp;
save -mat ~/Dropbox/temp/fitfunc.dat fit_X fit_Y
init_guess = [1e3 0 10e-6];
result = fminsearch('mygaussfit2', init_guess);
gauss_sigma_yp(n_bunch_qp) = result(3)
gauss_mu_yp = result(2)
gauss_A_yp = result(1)/gauss_sigma_yp(n_bunch_qp)/sqrt(2*pi);
xlabel('yp [rad]'); 



% energy
mean_E(n_bunch_qp) = mean(pp(:,6));
sigma_E_E(n_bunch_qp) = std(pp(:,6)) / mean_E(n_bunch_qp);


% calc correlation "cor" is for octave
CC = corrcoef( pp(:,5)*1e6, pp(:,1)*1e6);
zx_cor = CC(1,2)

% put in artifical x-z correlation to test
%zx_corr_art = 1 * 0
%zx_corr_art_2 = -10000* 0
%pp(:,1) = pp(:,1) + zx_corr_art*pp(:,5) + zx_corr_art_2*pp(:,5).^2;

% emittance calc
    Cxx=sqrt(det(cov(pp(:,[1,2]))));
    Cyy=sqrt(det(cov(pp(:,[3,4]))));
    emitt=[Cxx,Cyy]*mean(pp(:,6))/0.510998903076601*1e-6;
    emitt_x(n_bunch_qp) = emitt(1)
    emitt_y(n_bunch_qp) = emitt(2);

sigma_x(n_bunch_qp) = std(pp(:,1));
sigma_y(n_bunch_qp) = std(pp(:,3));
sigma_xp(n_bunch_qp) = std(pp(:,3));
sigma_yp(n_bunch_qp) = std(pp(:,4));
sigma_z(n_bunch_qp) = std(pp(:,5));


% emittance from gaussian fit - assuming no correlation .. (upright ellipse)
gauss_emitt_x = mean_E(n_bunch_qp)/0.511e6  .* gauss_sigma_x .* gauss_sigma_xp
gauss_emitt_y = mean_E(n_bunch_qp)/0.511e6  .* gauss_sigma_y .* gauss_sigma_yp


%
% find x, y correlation
%
tilt_x(:, n_bunch_qp) = polyfit(pp(:,5)*1e6, pp(:,1)*1e6,2)'
tilt_y(:, n_bunch_qp) = polyfit(pp(:,5)*1e6, pp(:,3)*1e6,2)'


if(0)

%min_Z = -626;
%max_Z = +626;

subplot(2,2,1);
my_densplot(my_ele2qp(pp), 3, 6, 0, hz_Nstd);
%plot( pp(:,5)*1e6, pp(:,6)/1e9, 'o');
%xlabel('z [um]'); 
%ylabel('E [GeV]');
%myaxis = axis;
%axis([min_Z  max_Z  myaxis(3)  myaxis(4)]);
%grid on;

subplot(2,2,3);
%[NN, XX] = hist(pp(:,5)*1e6, 51);
bar(XXz(n_bunch_qp, :)*1e6, NNz(n_bunch_qp, :));
xlabel('z [um]'); 
ylabel('count [a.u.]');
myaxis = axis;
axis([z_min*1e6  z_max*1e6  myaxis(3)  myaxis(4)]);
grid on;


subplot(2,2,2);
my_densplot(my_ele2qp(pp), 3, 1, 0, hz_Nstd);
%plot( pp(:,5)*1e6, pp(:,1)*1e6, 'o');
%xlabel('z [um]'); 
%ylabel('x [um]');
%myaxis = axis;
%axis([min_Z  max_Z  myaxis(3)  myaxis(4)]);
%axis([min_Z  max_Z  -600 400]);
grid on;
% together with estimated tilt
%hold on;
%plot(pp(:,5)*1e6, (pp(:,5)*1e6).^0*tilt_x(3) + (pp(:,5)*1e6).^1*tilt_x(2) + (pp(:,5)*1e6).^2*tilt_x(1), '.w');
%hold off;




subplot(2,2,4);
plot( pp(:,1)*1e6, pp(:,2)*1e6, 'o');
xlabel('x [um]'); 
ylabel('xp [urad]');
axis([-600 400 -50 40]);
axis([x_min  x_max  myaxis(3)  myaxis(4)]);
grid on;

subplot(2,2,4);
bar(XXx*1e6, NNx);
xlabel('x [um]'); 
ylabel('count [a.u.]');
grid on;

end% if

if( do_pause )
  pause;
end% if

end% for bunches

%
% OTHER INPUT TO QUICKPIC RPINPUT
%

% FACET parameters
plasma_s_prop = 1.0; % [m]  - plasma propagation length
Plasma_Density = 1.0e15; % /cm^3
%Plasma_Density = 2.0e16; % /cm^3
Plasma_Gas_Species = 3; % atom number
plasma_PREION=0; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step)
beam_match = 1; % override sigma_x, sigma_y with matched counterparts
emitt_match = 0;
%charge = charge * 2

% also, only input linear part
tilt_x([1 3], :) = 0;
tilt_y([1 3], :) = 0;

% force no tilt (for two bunches)
%tilt_x = [0 0 0; 0 0 0]';
tilt_y = [0 0 0; 0 0 0]';



%
% generate rpinput from elegant file
%
my_gen_rpinput(rpinput_template_file, rpinput_output_file, Plasma_Density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, gauss_sigma_x, gauss_sigma_y, sigma_z, gauss_emitt_x, gauss_emitt_y, BEAM_EV, beam_match, emitt_match, [1 1], tilt_x, tilt_y, XXz, NNz)




