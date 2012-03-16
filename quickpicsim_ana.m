%
% QuickPIC Matlab Analysis main script
% E. Adli, Dec 13, 2011
%
clear all;
clf;

%
%
%  user settings
%
%

% update working dir to your QUICKPICSIM folder
working_dir = '/Users/eadli/Dropbox/SLAC/quickpic/QUICKPICSIM/';

% update data dir to your quickpic output folder
datadir = '~/quicksimoutput/testrun/'



% choose output to analyze
do_QEB = 1;
do_QEB_3D = 0;
%do_QEB_line = 0;
do_QEP = 1;
%do_QEP_line = 0;
do_FEZ = 1;
%do_FEZ_line = 0;
% beam phase space analyze options
do_beam_phase_space = 1;
do_beam_phase_space_force_t0 = 0; % show phase space for first time step (only)
% the following assumes do_beam_phase_space = 1
do_beam_x_y = 0;
do_beam_x_xp = 0;
do_beam_y_yp = 0;
do_beam_z_E = 0;

% start plot  3D timestep, -1: start from start
%n_3D_start = 1800;
%n_3D_start = 600;
n_3D_start = -1;

 
% plotting?
do_plot = 1;
do_fixed_axes = 1;  % set manually below, if not auto-scale per graph
% for movie-makers
do_movie = 0;
n_movframe = 0;
N_dt_per_frame = 5; % for slowing down the movie, N frames for each pic


%
%
% end user settings
%
%


if(do_movie)
  set(0,'defaultaxesfontsize',13);
else
  set(0,'defaultaxesfontsize',18);
end% if

%
% general stuff
%
addpath(working_dir);
% import my standard SI constants
eval(['run ' working_dir 'my_SI_params.m']);

% color
mycolormap = 'jet';
c_min = 1e10;
c_max = -1e10;
% max for plotting
qp_QEB_max = -1e10;
qp_QEP_max = -1e10;
qp_FEZ_max = -1e10;
qp_QEB_min = 1e10;
qp_QEP_min = 1e10;
qp_FEZ_min = 1e10;


%
%
% import and convert quickpic data
%
%

% extract relevant data from RPINPUT file
myfile_rpinput = [datadir 'rpinput'];
n0 = my_get_quickpic_param(myfile_rpinput, 'Plasma_Density')
neutral_gas = my_get_quickpic_param(myfile_rpinput, 'Neutral_gas')
Box_X = my_get_quickpic_param(myfile_rpinput, 'Box_X') 
Box_Y = my_get_quickpic_param(myfile_rpinput, 'Box_Y') 
Box_Z = my_get_quickpic_param(myfile_rpinput, 'Box_Z') 
INDX = my_get_quickpic_param(myfile_rpinput, 'INDX') 
INDY = my_get_quickpic_param(myfile_rpinput, 'INDY') 
INDZ = my_get_quickpic_param(myfile_rpinput, 'INDZ') 
TEND = my_get_quickpic_param(myfile_rpinput, 'TEND')
DT = my_get_quickpic_param(myfile_rpinput, 'DT')
N_1 = my_get_quickpic_param(myfile_rpinput, 'Num_Particle') % for first beam
Gamma_1 = my_get_quickpic_param(myfile_rpinput, 'Gamma') % for first beam
DFQEBSLICE = my_get_quickpic_param(myfile_rpinput, 'DFQEBSLICE')
DFQEPSLICE = my_get_quickpic_param(myfile_rpinput, 'DFQEBSLICE')
DFESLICE = my_get_quickpic_param(myfile_rpinput, 'DFESLICE')
if( n0 == -1 | Box_X == -1 | Box_Y == -1 | Box_Z == --1 | INDX == -1 | INDY == -1 | INDZ == -1 | TEND == -1 | DT == -1 | DFQEBSLICE == -1 | DFQEPSLICE == -1 | DFESLICE == -1)
  warning('EA: a rpinput parameter could not be read.  Stopping execution.');
  stop;
end% if
if( (DFQEBSLICE ~= DFQEPSLICE) | (DFQEBSLICE ~= DFESLICE) | (DFESLICE ~= DFQEPSLICE) )
  warning('EA: slice interval set differently for different outout.  Using DFQEBSLICE for timing.');
end% if


% LOOP 3D time step
n_3D_counter = 0;
N_3D_timestep = floor(TEND/DT)
if(n_3D_start < 0)
  n_3D_timestep_start = DFQEBSLICE;
else
  n_3D_timestep_start = n_3D_start;
end% if
for n_3D_timestep = n_3D_timestep_start:DFQEBSLICE:N_3D_timestep,
n_3D_timestep_str = sprintf('%.4d', n_3D_timestep');
n_3D_counter = n_3D_counter + 1;
% skip update for beam to store all data, took too much memory
%n_3D_counter_beam = n_3D_counter;
n_3D_counter_beam = 1;

%
% read specified field and dists
%
% establish quickpic output format
qp_version_suffix = my_get_quickpic_format(datadir);
if (length (qp_version_suffix) == 0)
  disp('EA: unknown QuickPIC version format');
  stop;
end% if
if(do_FEZ)
  myfile = [datadir 'FEZ-XZ/FEZ-XZ_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).FEZ = double(my_read_hdf(myfile));
end% if
if(do_QEB)
  myfile = [datadir 'QEB-XZ/QEB-XZ_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).QEB = double(my_read_hdf(myfile));
end% if
if(do_QEB_3D)
  myfile = [datadir 'QEB/QEB_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).QEB_3D = double(my_read_hdf(myfile));
end% if
if(do_QEP)
  if( strfind(qp_version_suffix, '.h5') )
    myfile = [datadir 'QEP1-XZ/QEP1-XZ_' n_3D_timestep_str qp_version_suffix];
  else
    myfile = [datadir 'QEP01-XZ/QEP01-XZ_' n_3D_timestep_str qp_version_suffix];
  end% if
  qp(n_3D_counter).QEP = double(my_read_hdf(myfile));
end% if

%
% scale data to physical units
%
scale_x = Box_X / 2^INDX;
scale_y = Box_Y / 2^INDY;
scale_z = Box_Z / 2^INDZ;
offset_x0 = Box_X/2; % assume initial beam is put in the middle of the box
offset_y0 = Box_Y/2; % assume initial beam is put in the middle of the box
offset_z0 = Box_Z/2; % assume initial beam is put in the middle of the box
scale_E = 1;
scale_rho = 1;
% norm to cgs
omega_p = sqrt(n0*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
scale_E = scale_E * SI_em*1e3*SI_c*1e2*omega_p / (SI_e / 3.336e-10);
scale_rho = scale_rho * (SI_e / 3.336e-10) * n0*1e6;
% cgs to mks
scale_E = scale_E / ( 1e4 / SI_c);
scale_rho = scale_rho / (1e-5*SI_c);
% V/m to GV/m
scale_E = scale_E / 1e9;
% C/m^3 to N / cm^3
scale_rho = scale_rho/SI_e / 1e6;
% N / cm^3 to units of [np_0]
scale_rho = scale_rho / n0/1e6;
% update data
if(do_FEZ)
  qp(n_3D_counter).FEZ = qp(n_3D_counter).FEZ * scale_E;
end% if
if(do_QEB)
  qp(n_3D_counter).QEB = qp(n_3D_counter).QEB * scale_rho * -1;
end% if
if(do_QEB_3D)
  qp(n_3D_counter).QEB_3D = qp(n_3D_counter).QEB_3D * scale_rho * -1;
  % store total beam charge
  qp(n_3D_counter).Q_total = sum(sum(sum(qp(n_3D_counter).QEB_3D)));
end% if
if(do_QEP)
  qp(n_3D_counter).QEP = qp(n_3D_counter).QEP * scale_rho * -1;
end% if


% read beam phase space
if(do_beam_phase_space )

n_beams = my_get_quickpic_phasespace_nbeams(datadir); % total # of beams
for n_beam=1:n_beams,
% analyze each beam
n_beam_str = sprintf('%.2d', n_beam');
if(do_beam_phase_space_force_t0 )
  n_3D_timestep = 0; % override, to look at a single time step dump
end% if
%  n_3D_timestep = 1240; % override, to look at a single time step dump
n_3D_timestep_str = sprintf('%.4d', n_3D_timestep');

% load phase space beam (merge all beam parts)
qp(n_3D_counter_beam).PP(n_beam).BEAM = my_get_quickpic_phasespace(datadir, n_beam_str, n_3D_timestep_str); % auto-extracts beam part filenames exists
% set units
qp_version_suffix = my_get_quickpic_format(datadir);
if( strfind(qp_version_suffix, '.h5') )
  % scaling for newer quickpic versions
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,1) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,1) * 1 / k_p * 1e6; % to um
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,2) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,2) * 1 / k_p * 1e6; % to um
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,3) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,3) * 1 / k_p * 1e6; % to um
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,4) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,4) ./ ( mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6)) ) * 1e6; % to urad
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,5) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,5) ./ ( mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6)) ) * 1e6; % to urad
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6) * SI_em * SI_c^2/SI_e / 1e9; % total energy [GeV]
else
  % scaling for older quickpic versions
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,1) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,1) * scale_x - offset_x0; % to um
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,2) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,2) * scale_y - offset_y0; % to um
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,3) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,3) * scale_z - offset_z0; % to um
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,4) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,4) ./ ( mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6)) ) * 1e6; % to urad
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,5) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,5) ./ ( mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6)) ) * 1e6; % to urad
  qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6) = qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6) * SI_em * SI_c^2/SI_e / 1e9; % total energy [GeV]
end% 

% store some quick numbers
qp(n_3D_counter).PP(n_beam).sigma_x = std(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,1));
qp(n_3D_counter).PP(n_beam).sigma_y = std(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,2));
qp(n_3D_counter).PP(n_beam).sigma_z = std(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,3));
qp(n_3D_counter).PP(n_beam).sigma_xp = std(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,4));
qp(n_3D_counter).PP(n_beam).sigma_yp = std(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,5));
qp(n_3D_counter).PP(n_beam).std_E = std(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6));
qp(n_3D_counter).PP(n_beam).mean_x = mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,1));
qp(n_3D_counter).PP(n_beam).mean_y = mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,2));
qp(n_3D_counter).PP(n_beam).mean_z = mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,3));
qp(n_3D_counter).PP(n_beam).mean_xp = mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,4));
qp(n_3D_counter).PP(n_beam).mean_yp = mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,5));
qp(n_3D_counter).PP(n_beam).mean_E = mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6));
qp(n_3D_counter).PP(n_beam).emitt = [sqrt(det(cov(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,[1,4]) * 1e-6 ))), sqrt(det(cov(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,[2,5])* 1e-6 )))] * mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6)) / ( SI_em * SI_c^2/SI_e) *1e9;
qp(n_3D_counter).FEZ_max = max(max(qp(n_3D_counter).FEZ(1:round(end*3/3), 1:end)));;
qp(n_3D_counter).FEZ_min = min(min(qp(n_3D_counter).FEZ(1:round(end*3/3), 1:end)));
qp(n_3D_counter).QEB_max =  max(max(qp(n_3D_counter).QEB));
qp(n_3D_counter).PP(n_beam).mean_E = mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6));
% store slice transverse data
if(n_3D_counter == 1)
  [slice_mean_x, slice_sigma_x, slice_z, slice_N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 1);
else
  % ensure same slicing as beam evolves through plasma, in order to compare apples to apples
  [slice_mean_x, slice_sigma_x, slice_z, slice_N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 1, slice_z);
end% if
qp(n_3D_counter).PP(n_beam).slice_z = slice_z;
qp(n_3D_counter).PP(n_beam).slice_N_z = slice_N_z;
qp(n_3D_counter).PP(n_beam).slice_mean_x = slice_mean_x;
qp(n_3D_counter).PP(n_beam).slice_sigma_x = slice_sigma_x;
% store slice energy data
if(n_3D_counter == 1)
  [slice_mean, slice_sigma, slice_z, slice_N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 6);
else
  [slice_mean, slice_sigma, slice_z, slice_N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 6, slice_z);
end% if
qp(n_3D_counter).PP(n_beam).slice_z = slice_z;
qp(n_3D_counter).PP(n_beam).slice_N_z = slice_N_z;
qp(n_3D_counter).PP(n_beam).slice_mean_E = slice_mean;
qp(n_3D_counter).PP(n_beam).slice_sigma_E = slice_sigma;
% store slice twiss data
[emnx,emny,betax,alphax,betay,alphay] = my_get_slice_twiss(qp(n_3D_counter_beam).PP(n_beam).BEAM, slice_z);
qp(n_3D_counter).PP(n_beam).slice_emnx = emnx;
qp(n_3D_counter).PP(n_beam).slice_emny = emny;
qp(n_3D_counter).PP(n_beam).slice_betax = betax;
qp(n_3D_counter).PP(n_beam).slice_alphax = alphax;
qp(n_3D_counter).PP(n_beam).slice_betay = betay;
qp(n_3D_counter).PP(n_beam).slice_alphay = alphay;
end% for each beam

qp(n_3D_counter).n0 = n0;
s_timestep = (SI_c/omega_p) * DT * n_3D_timestep; % propagation length for this timestep
qp(n_3D_counter).s_timestep = s_timestep;
disp(' ' ); 
display(['Prop. length s=' num2str(s_timestep*100, 3) ' [cm]. Time steps in units of DT: ' n_3D_timestep_str]);

end% if phase space


         % put in this to separate data read and plot; however,
         % data size becomes to large
         %end% 3D timestep loop
% LOOP 3D time step
%n_3D_counter = 0;
%for n_3D_timestep = DFQEBSLICE:DFQEBSLICE:N_3D_timestep,
%for n_3D_timestep = 2220:DFQEBSLICE:N_3D_timestep,
%n_3D_timestep_str = sprintf('%.4d', n_3D_timestep');
%n_3D_counter = n_3D_counter + 1;

         
         
%         
%         
% calc some plasma quantities
%
%

if( do_beam_phase_space )
I_peak = N_1*SI_e *SI_c/sqrt(2*pi*(qp(n_3D_counter).PP(n_beam).sigma_z/1e6)^2);
n_dens_bunch = N_1 / ( (2*pi)^(3/2)*qp(n_3D_counter).PP(n_beam).sigma_x/1e6*qp(n_3D_counter).PP(n_beam).sigma_y/1e6*qp(n_3D_counter).PP(n_beam).sigma_z/1e6 ) / 1e6; % cm^-3
f_betatron = sqrt(2*Gamma_1)*(2*pi) / k_p;
Lambda_nonlin = 2.5*(N_1/2e10)*(20/qp(n_3D_counter).PP(n_beam).sigma_x);
crit_blow_out = n_dens_bunch / (n0 * (1+4/(k_p*qp(n_3D_counter).PP(n_beam).sigma_z/1e6)^2) );
V_erosion_PI = sqrt( qp(n_3D_counter).PP(n_beam).emitt(1)*k_p/(qp(n_3D_counter).PP(n_beam).mean_E/.511e-3)^(3/2)/(I_peak/17e3) );
E0_wb = SI_em*SI_c*omega_p / SI_e / 1e9; % Wave-breaking field, GeV/m
Ez_lin_IB = 4*pi*E0_wb * I_peak / 17e3;
Ez_lin = -100*(N_1/2e10)*(20/qp(n_3D_counter).PP(n_beam).sigma_z)^2*log(sqrt(2.5e17*1e6 / n0/1e6 * 10 / qp(n_3D_counter).PP(n_beam).sigma_x)); % lin. regime (NJP), GeV/m
Ez_sel = 40*(N_1/2e10)*(15/qp(n_3D_counter).PP(n_beam).sigma_x)*(20/qp(n_3D_counter).PP(n_beam).sigma_z)^2;
Ez_hog = 0.244*(N_1/2e10)*(600/qp(n_3D_counter).PP(n_beam).sigma_z)^2;
end% if


%
%
% plot data
%
%

if(do_plot)

clf;

%
%  plot QEB
%
if(do_QEB)
  subplot(2,4,2);
  colormap(mycolormap);
  ZZ = (1:size(qp(n_3D_counter).QEB,1)) * scale_z - offset_z0;
  XX = (1:size(qp(n_3D_counter).QEB,2)) * scale_x - offset_x0;
  YY = (1:size(qp(n_3D_counter).QEB,2)) * scale_y - offset_y0;
  c_min_new = (min( min(min(qp(n_3D_counter).QEB)), min(min(qp(n_3D_counter).QEB)) ));
  c_max_new = (max( max(max(qp(n_3D_counter).QEB)), max(max(qp(n_3D_counter).QEB)) ));
  % don't jump up and down as time evolve
  if(c_min_new < c_min)
    c_min = c_min_new;
  end% if
  if(c_max_new > c_max)
    c_max = c_max_new;
  end% if

%  pcolor(ZZ,XX,qp(n_3D_counter).QEB');
  pcolor(ZZ,XX,qp(n_3D_counter).QEB'+qp(n_3D_counter).QEP');  % with plasma
%  pcolor(ZZ,XX,qp(n_3D_counter).QEB');  % with plasma
  h_c = colorbar('NorthOutside');
  %caxis([0.0 0.1]);
  if(do_fixed_axes)
    caxis([0 100]); % NJP
    caxis([0 30]);
    caxis([0 5]); % Dx-3
    caxis([0 5]); % Dx-0 mat
    caxis([-10 10]); % positrons
  end% if
    caxis([0 10]); % positrons
  shading('flat');
  xlabel('z [um]')
  ylabel('x [um]')
  %set(get(h_c,'ylabel'),'String', 'n_b  [n_p]');
  %set(get(h_c,'ylabel'),'String', 'n_b  [C/cm^3]');
%  set(get(h_c,'ylabel'),'String', 'n_b  [n_p_0]');
  set(get(h_c,'ylabel'),'String', 'n_e  [n_p_0]');  % with plasma
  
  % add line plot on axis
  subplot(2,4,6);
  plot(ZZ, qp(n_3D_counter).QEB(:,(2^INDX/2)+1)); % middle line
%  plot(ZZ, qp(n_3D_counter).QEB(:,(2^INDX/2)-5)); % middle line
  %plot(ZZ, sum(qp(n_3D_counter).QEB, 2)); % integrated
  grid on;
  xlabel('z [um]')
  ylabel('n_b [n_p_0]')
  qp_QEB_max = max(qp_QEB_max, max(qp(n_3D_counter).QEB(:,(2^INDX/2)+1))); % middle line 
  qp_QEB_min = min(qp_QEB_min, min(qp(n_3D_counter).QEB(:,(2^INDX/2)+1))); % middle line
  %qp_QEB_max = max(qp_QEB_max, max( sum(qp(n_3D_counter).QEB,2))); % integrated
  %qp_QEB_min = min(qp_QEB_min, min( sum(qp(n_3D_counter).QEB,2))); % integrated
  myaxis = axis;
  axis([min(ZZ)  max(ZZ)  qp_QEB_min  qp_QEB_max+eps]);
  QEB_max_0 =  max(max(qp(1).QEB));
  if(do_fixed_axes)
    axis([min(ZZ)  max(ZZ)  0 150]); % NJP
    axis([min(ZZ)  max(ZZ)  0 100]);
    %axis([min(ZZ)  max(ZZ)  0 QEB_max_0]); % Dx03
    axis([min(ZZ)  max(ZZ)  -40 40]); % positrons
  end% if
  QEB_max =  max(max(qp(n_3D_counter).QEB));
  title(['n_{b,max} = ' num2str(QEB_max, '%.1f') ' n_p_0']);
%  pause;
end% if plot


%
%  plot QEP
%
if(do_QEP)
  subplot(2,4,3);
  colormap(mycolormap);
  ZZ = (1:size(qp(n_3D_counter).QEP,1)) * scale_z - offset_z0;
  XX = (1:size(qp(n_3D_counter).QEP,2)) * scale_x - offset_x0;
  YY = (1:size(qp(n_3D_counter).QEP,2)) * scale_y - offset_y0;
  c_min_new = (min( min(min(qp(n_3D_counter).QEP)), min(min(qp(n_3D_counter).QEP)) ));
  c_max_new = (max( max(max(qp(n_3D_counter).QEP)), max(max(qp(n_3D_counter).QEP)) ));
  % don't jump up and down as time evolve
  if(c_min_new < c_min)
    c_min = c_min_new;
  end% if
  if(c_max_new > c_max)
    c_max = c_max_new;
  end% if

  pcolor(ZZ,XX,qp(n_3D_counter).QEP');
  h_c = colorbar('NorthOutside');
%  caxis([c_min c_max]);
  caxis([0 5]);
  if(do_fixed_axes)
    caxis([0 2]);
    caxis([0 10]); % positrons
  end% if
  shading('flat');
  xlabel('z [um]')
  ylabel('x [um]')
  %set(get(h_c,'ylabel'),'String', 'n_p  [n_p]');
  %set(get(h_c,'ylabel'),'String', 'n_p  [C/cm^3]');
  set(get(h_c,'ylabel'),'String', 'n_p  [n_p_0]');

  % add line plot on axis
  subplot(2,4,7);
  plot(ZZ, qp(n_3D_counter).QEP(:,(2^INDX/2)+1));
  grid on;
  xlabel('z [um]')
  ylabel('n_p [n_p_0]')
  qp_QEP_max = max(qp_QEP_max, max(qp(n_3D_counter).QEP(:,(2^INDX/2)+1)));
  qp_QEP_min = min(qp_QEP_min, min(qp(n_3D_counter).QEP(:,(2^INDX/2)+1)));
  myaxis = axis;
  axis([min(ZZ)  max(ZZ)  qp_QEP_min  qp_QEP_max+eps]);
  axis([min(ZZ)  max(ZZ)  qp_QEP_min 2]);
  if(do_fixed_axes)
    axis([min(ZZ)  max(ZZ)  0 2]);
    axis([min(ZZ)  max(ZZ)  0 10]); % positrons
  end% if
  title(['n_p_0=' num2str(n0, '%.1e') '/cm^3' '[' num2str(neutral_gas) ']']);

%  pause;
end% if plot



%
%  plot FEZ
%
if(do_FEZ)
  subplot(2,4,4);
  colormap(mycolormap);
  ZZ = (1:size(qp(n_3D_counter).FEZ,1)) * scale_z - offset_z0;
  XX = (1:size(qp(n_3D_counter).FEZ,2)) * scale_x - offset_x0;
  YY = (1:size(qp(n_3D_counter).FEZ,2)) * scale_y - offset_y0;
  c_min_new = (min( min(min(qp(n_3D_counter).FEZ)), min(min(qp(n_3D_counter).FEZ)) ));
  c_max_new = (max( max(max(qp(n_3D_counter).FEZ)), max(max(qp(n_3D_counter).FEZ)) ));
  % don't jump up and down as time evolve
  if(c_min_new < c_min)
    c_min = c_min_new;
  end% if
  if(c_max_new > c_max)
    c_max = c_max_new;
  end% if

  pcolor(ZZ,XX,qp(n_3D_counter).FEZ');
  h_c = colorbar('NorthOutside');
%  caxis([c_min c_max]);
  if(do_fixed_axes)
    caxis([-70 70]); % NJP
    caxis([-30 30]);
    caxis([-5 5]); % Dx-3
    caxis([-5 5]); % positrons
  end% if
  shading('flat');
  xlabel('z [um]')
  ylabel('x [um]')
  %set(get(h_c,'ylabel'),'String', 'n_b  [n_p]');
  %set(get(h_c,'ylabel'),'String', 'n_b  [C/cm^3]');
  set(get(h_c,'ylabel'),'String', 'E_z [GV/m]');

  % add line plot on axis
  subplot(2,4,8);
  plot(ZZ, qp(n_3D_counter).FEZ(:,(2^INDX/2)+1));
  grid on;
  xlabel('z [um]')
  ylabel('E_z [GV/m]')
  qp_FEZ_max = max(qp_FEZ_max, max(qp(n_3D_counter).FEZ(:,(2^INDX/2)+1)));
  qp_FEZ_min = min(qp_FEZ_min, min(qp(n_3D_counter).FEZ(:,(2^INDX/2)+1)));
  myaxis = axis;
  axis([min(ZZ)  max(ZZ)  qp_FEZ_min  qp_FEZ_max+eps]);
  FEZ_max_0 = max(max(qp(1).FEZ(1:round(end*2/3), 1:end)));
  if(do_fixed_axes)
    axis([min(ZZ)  max(ZZ) -70 70]); % NJP
    axis([min(ZZ)  max(ZZ) -50 50]);
    axis([min(ZZ)  max(ZZ) -FEZ_max_0 FEZ_max_0]); % Dx-3
  end% if
  FEZ_max = max(max(qp(n_3D_counter).FEZ(1:round(end*2/3), 1:end)));
  FEZ_min = min(min(qp(n_3D_counter).FEZ(1:round(end*2/3), 1:end)));
  title(['E_{z,max} = ' num2str(FEZ_max, '%.1f') 'GeV/m']);
  %title(['E_{z,max} = ' num2str(FEZ_max, '%.1f') 'GeV/m'  ', ' 'E_{z,min} = ' num2str(FEZ_min, '%.1f') 'GeV/m']);

%  pause;
end% if plot



%
%  plot beam phase space
%
if(do_beam_phase_space )

% total beam analysis - merge beams (not stored in main qp var for space reasons)
qp_BEAMS = [];
for(n_beam=1:size(qp(n_3D_counter_beam).PP, 2))
  qp_BEAMS = [qp_BEAMS; qp(n_3D_counter_beam).PP(n_beam).BEAM];
end% for

% x, E histogram;  x : put 5 initial sigma
subplot(2,4,[1 5]);
%hist_var1 = 3; % abscissa
hist_var1 = 1; % abscissa
hist_var2 = 4; % ordinate
pplabel(1).var = 'x [um]';
pplabel(2).var = 'y [um]';
pplabel(3).var = 'z [um]';
pplabel(4).var = 'xp [urad]';
pplabel(5).var = 'yp [urad]';
pplabel(6).var = 'p [GeV/c]';
%x_min = min(qp_BEAMS(:,hist_var1));
%x_max = max(qp_BEAMS(:,hist_var1));
x_min = mean(qp_BEAMS(:,hist_var1))-5*std(qp_BEAMS(:,hist_var1));
x_max = mean(qp_BEAMS(:,hist_var1))+5*std(qp_BEAMS(:,hist_var1));
x_min = mean(qp_BEAMS(:,hist_var1))-5*std(qp_BEAMS(:,hist_var1)); % z-E plot
x_max = mean(qp_BEAMS(:,hist_var1))+5*std(qp_BEAMS(:,hist_var1)); % z-E plot
if(do_fixed_axes)
  x_min = -40;
  x_max = +40;
  if(hist_var1 == 3)
    x_min = -0;  % NJP
    x_max = +300;% NJP
    x_min = -330;
    x_max = +0;
    x_min = -330; % Dx-3
    x_max = +330;   % Dx-4
  end    
end% if
E_min = min(qp_BEAMS(:,hist_var2));
E_max = max(qp_BEAMS(:,hist_var2));
E_min
if(do_fixed_axes)
  E_min = 0;
  E_max = 60;
  E_min = -40;
  E_max = 40;
  E_min = 15; % Dx-3
  E_max = 30;% Dx-3
  E_min = 40; % positrons
  E_max = 60;% positrons
  E_min = 0; % NJP
  E_max = 60;% NJP
end% if
x_min = -50;
x_max= 50;
  E_min = -10000; % x xp
  E_max = 10000;% x xp
xedges = linspace(x_min, x_max, 512+0*round((x_max-x_min)*scale_x));
yedges = linspace(E_min, E_max, 512);
histmat = hist2(qp_BEAMS(:,hist_var1), qp_BEAMS(:,hist_var2), xedges, yedges);
h_g = pcolor(xedges,yedges,histmat); 
h_c = colorbar;
if(do_fixed_axes)
  caxis([1 100]); % NJP
  caxis([0 10]);
  caxis([0 20]); % positrons
end% if
%  caxis([0 10]);
shading('flat');
axis([x_min x_max E_min E_max]); % for comp with AAC 2010 figure
%axis square tight;
xlabel(pplabel(hist_var1).var);
ylabel(pplabel(hist_var2).var);
%set(get(h_c,'ylabel'),'String', 'n_b  [a.u.]');
title(['s=' num2str(s_timestep*100, 3) ' [cm].   Time step: ' num2str(n_3D_timestep) ' [DT]']);

%pause;

% movie-maker
if( do_movie )
for(n_dt_dt_per_frame=1:N_dt_per_frame)
  n_movframe = n_movframe + 1;
  mymovdir = [working_dir 'movie_frames/'];
  myfile_mov = [mymovdir 'frame' num2str(n_movframe, '%5.5d') '.png'];
  saveas(h_g, myfile_mov, 'png');
end% for
end% if

%  display(['Beam ' n_beam_str ': ' 'sigma_x=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_x(n_3D_counter_beam)) '  sigma_y=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_y(n_3D_counter_beam)) '  sigma_z=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_z(n_3D_counter_beam))  '  sigma_xp=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_xp(n_3D_counter_beam))  '  sigma_yp=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_yp(n_3D_counter_beam))   '  mean_E=' num2str(qp(n_3D_counter_beam).PP(n_beam).mean_E(n_3D_counter_beam))     '  std_E / mean_E=' num2str(qp(n_3D_counter_beam).PP(n_beam).std_E(n_3D_counter_beam) / qp(n_3D_counter_beam).PP(n_beam).mean_E(n_3D_counter_beam)) ]);

if(~do_movie)
  pause;
end% if

end% if plot
    
end% 3D timestep loop

end% if plot


