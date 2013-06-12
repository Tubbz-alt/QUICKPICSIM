%
% QuickPIC Matlab Analysis main script
% E. Adli, Dec 13, 2011
%
% Last update: E. Adli, Jun 12, 2013
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

%
% CONVERGENCE TEST
datadir = '/Users/eadli/templ/ctm/FI_sr1_3/';
datadir = '/Users/eadli/templ/ctm/FI_sr1_3_999/';
datadir = '/Users/eadli/templ/ctm/FI_sr1_3_101010/';
% plus/minus in box scan
datadir = '/Users/eadli/templ/ctm/sr4/';
datadir = '/Users/eadli/templ/ctm/sr4_p5/';
datadir = '/Users/eadli/templ/ctm/sr4_m5/';
% TILT SCAN (FI)
datadir = '/Users/eadli/templ/ctm/sr4_t000/';
datadir = '/Users/eadli/templ/ctm/sr4/';
datadir = '/Users/eadli/templ/ctm/sr4_t010/';
datadir = '/Users/eadli/templ/ctm/sr4_t020/';
%PI_sr1_3___t0025 - t010 - tilt tests
datadir = '/Users/eadli/templ/ctm/PI_sr1_3___t0025/';
datadir = '/Users/eadli/templ/ctm/PI_sr1_3___t010/';
datadir = '/Users/eadli/templ/ctm/PI_sr1_3/';
%PI_sr1_3___t001 - convergence tests
datadir = '/Users/eadli/templ/ctm/PI_sr1_3___t001/';
%datadir = '/Users/eadli/templ/ctm/PI_sr1_3___t001___999/';
%datadir = '/Users/eadli/templ/ctm/PI_sr1_3___t001___101010/';
% EMITT SCAN unmatched, const size - FI
datadir = '/Users/eadli/templ/ctm/FI_sr5_em5/';
%datadir = '/Users/eadli/templ/ctm/FI_sr5_em4/';
%datadir = '/Users/eadli/templ/ctm/FI_sr5_em3/';
%datadir = '/Users/eadli/templ/ctm/FI_sr5_em1_8/';
%datadir = '/Users/eadli/templ/ctm/FI_sr5_em1_3/';
%datadir = '/Users/eadli/templ/ctm/FI_sr5_em0_9/';
%datadir = '/Users/eadli/templ/ctm/FI_sr5_em0_1/';
%datadir = '/Users/eadli/templ/ctm/FI_sr5_em0_01/';
% EMITT SCAN unmatched, const size - PI
datadir = '/Users/eadli/templ/ctm/PI_sr5_em5/';
%datadir = '/Users/eadli/templ/ctm/PI_sr5_em4/';
%datadir = '/Users/eadli/templ/ctm/PI_sr5_em3/';
%datadir = '/Users/eadli/templ/ctm/PI_sr5_em1_8/';
%datadir = '/Users/eadli/templ/ctm/PI_sr5_em1_3/';
%datadir = '/Users/eadli/templ/ctm/PI_sr5_em0_9/';
%datadir = '/Users/eadli/templ/ctm/PI_sr5_em0_1/';
%datadir = '/Users/eadli/templ/ctm/PI_sr5_em0_01/';
% EMITT SCAN (matched) - PI
%datadir = '/Users/eadli/templ/ctm/PI_sr5/';
%datadir = '/Users/eadli/templ/ctm/PI_sr4/';
%datadir = '/Users/eadli/templ/ctm/PI_sr3/';
%datadir = '/Users/eadli/templ/ctm/PI_sr1_8/';
%datadir = '/Users/eadli/templ/ctm/PI_sr1_3/';
%datadir = '/Users/eadli/templ/ctm/PI_sr0_9/';
%datadir = '/Users/eadli/templ/ctm/PI_sr0_3/';
%datadir = '/Users/eadli/templ/ctm/PI_sr0_1/';
%datadir = '/Users/eadli/templ/ctm/PI_sr0_04/';
datadir = '/Users/eadli/templ/ctm/PI_sr0_01/';
% EMITT SCAN (matched) - FI
%datadir = '/Users/eadli/templ/ctm/FI_sr5/';
%datadir = '/Users/eadli/templ/ctm/FI_sr4/';
%datadir = '/Users/eadli/templ/ctm/FI_sr3/';
%datadir = '/Users/eadli/templ/ctm/FI_sr1_8/';
%datadir = '/Users/eadli/templ/ctm/FI_sr1_3/';
%datadir = '/Users/eadli/templ/ctm/FI_sr0_9/'
%datadir = '/Users/eadli/templ/ctm/FI_sr0_3/';
%datadir = '/Users/eadli/templ/ctm/FI_sr0_1/';
%datadir = '/Users/eadli/templ/ctm/FI_sr0_04/';
%datadir = '/Users/eadli/templ/ctm/FI_sr0_01/';

datadir = '/Users/eadli/templ/ctm/PI_sr1_3/';
datadir = '/Users/eadli/templ/ctm/PI_sr1_3___t001/';
datadir = '/Users/eadli/templ/ctm/PI_sr1_3___t001___999/';
datadir = '/Users/eadli/templ/ctm/FI_sr4/';

% LC OPTIM: T=1.0
datadir = '/Users/eadli/templ/LC/O2wb/'
datadir = '/Users/eadli/templ/LC/O1long/'
datadir = '/Users/eadli/templ/LC/O2long/'

% Ken SB
datadir = '/Users/eadli/templ/SB/1e17/'
datadir = '/Users/eadli/templ/SB/2e17/'
datadir = '/Users/eadli/templ/SB/4e17/'
datadir = '/Users/eadli/templ/SB/2e17_40/'
datadir = '/Users/eadli/templ/SB/2e17_60/'

% LC OPTIM: T=1.0
datadir = '/Users/eadli/templ/LC/O2wb/'
%datadir = '/Users/eadli/templ/LC/O2long/'

% exp hosing 2013
% beta 5 x 0.5 m2
datadir = '/Users/eadli/templ/hose_exp/beta5/2e17_t0/'
datadir = '/Users/eadli/templ/hose_exp/beta5/2e17_t05/'
datadir = '/Users/eadli/templ/hose_exp/beta5/4e17_t05_em2/'
datadir = '/Users/eadli/templ/hose_exp/beta5/5e16_t05_em2/'
datadir = '/Users/eadli/templ/hose_exp/beta5/5e16_t05/'
datadir = '/Users/eadli/templ/hose_exp/beta5/2_7e17_40_t05/'
datadir = '/Users/eadli/templ/hose_exp/beta5/1_1e17_60_t05/'
datadir = '/Users/eadli/templ/hose_exp/beta5/6_6e16_60_t05/'
% beta 0.5 x 0.06 m2
datadir = '/Users/eadli/templ/hose_exp/beta05/2_7e17_40_t05/'
% beta 1 x 0.1 m2
datadir = '/Users/eadli/templ/hose_exp/beta1/2_7e17_40_t05/'
datadir = '/Users/eadli/templ/hose_exp/beta1/1_1e17_40_t05/'
datadir = '/Users/eadli/templ/hose_exp/beta1/6_6e16_60_t05/'
datadir = '/Users/eadli/templ/hose_exp/beta1/1_1e17_60_t05/'
datadir = '/Users/eadli/templ/hose_exp/beta1/1_1e17_40_t05_10/'
%datadir = '/Users/eadli/templ/hose_exp/beta1/1_1e17_60_t05_10/'
datadir = '/Users/eadli/templ/hose_exp/beta1/2_7e17_40_t05_10/'
datadir = '/Users/eadli/templ/hose_exp/beta1/2_7e17_40_t01_10/'
datadir = '/Users/eadli/templ/hose_exp/beta1/1_1e17_60_t05_10/'
datadir = '/Users/eadli/templ/hose_exp/beta1/2_7e17_60_t01_10/'


datadir = '/Users/eadli/templ/ctm/FI_sr5/';
datadir = '~/templ/oldqp/n3e17_tilt001/';
datadir = '~/templ/dumpwork/3e17_se/';


%
%
% >>> USER SETTINGS
%
%


% plotting?
do_plot = 1;


% start plot  3D timestep, -1: start from start
n_3D_start = -1;
%n_3D_start = 1400;

% -1: go all the way to end
n_3D_end = -1;
%n_3D_end = 300;

% choose output to analyze
do_QEB_3D = 0;
do_FALL = 1;
% beam phase space analyze options
do_beam_phase_space = 1;
do_beam_phase_space_force_t0 = 0; % show phase space for first time step (only)
do_ana_slice__gauss_fit = 0; % [for all slices] time consuming, so can be taken out

% choose output to plot
do_FULL_beam_phase_space_plot = 0;
do_SLICE_beam_phase_space_plot = 0;
do_both_plot_types = 1; % only if above is 1
%n_slice_plot = 21; % z=-20
%n_slice_plot = 31; % z=20
n_slice_plot = 26; % z=0
do_calc_hose_ctm = 1;
do_QEB = 1;
do_QEP = 1;
do_FEZ = 1; 
do_FEX = 0;
do_FBY = 0;
do_FEY = 0;
do_FBX = 0;
do_FBZ = 0;
do_Eprof = 0;
do_Eprof_allview = 0;
do_EcrossB = 0;
do_FOCUSAXIS = 0;
do_QEB_FEZ = 0;
do_beam_phase_space_plot = 1;
do_translopeEB = 1;
translopeEB_dz = -3; % dz [um] from center of box
%translopeEB_dz = 50; % dz [um] from center of box
%translopeEB_dz = 'ion_head'; % if 'ion_head' use (z1_ion_head + z_ion_head)/2 to calculate z-slice of EB calc
%translopeEB_dz = 'hose_tail'; % if 'hose_tail' use (z1_tail_start + z_tail_end)/2 to calculate z-slice of EB calc

% more plot options
do_fixed_axes = 0;  % set manually below, if not auto-scale per graph
% for movie-makers
do_movie = 0;
n_movframe = 0;
N_dt_per_frame = 5; % for slowing down the movie, N frames for each pic


%
%
% <<< END USER SETTINGS
%
%



set(gcf, 'Color', 'w');
if(do_movie)
  set(0,'defaultaxesfontsize',18);
else
  set(0,'defaultaxesfontsize',24);
end% if

%  set(0,'defaultaxesfontsize',12);

colorlist = [{'r', 'b', 'g', 'm', 'c', 'r', 'b', 'g', 'm', 'c'}];

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
qp_FEX_max = -1e10;
qp_FBY_max = -1e10;
qp_FEY_max = -1e10;
qp_FBX_max = -1e10;
qp_FBZ_max = -1e10;
qp_QEB_min = 1e10;
qp_QEP_min = 1e10;
qp_FEZ_min = 1e10;
qp_FEX_min = 1e10;
qp_FBY_min = 1e10;
qp_FEY_min = 1e10;
qp_FBX_min = 1e10;
qp_FBZ_min = 1e10;
Eprof_n_max_glob = 0;


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
NBeams = my_get_quickpic_param(myfile_rpinput, 'NBeams')
for n = 1:NBeams,
  qp_N(n) = my_get_quickpic_param(myfile_rpinput, 'Num_Particle', n)
  qp_Gamma(n) = my_get_quickpic_param(myfile_rpinput, 'Gamma', n) 
  X0(n) = my_get_quickpic_param(myfile_rpinput, 'Parameter_Array(1:1,1:3)', n);
end% for
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
if(n_3D_start < 0)
  n_3D_timestep_start = 0;
else
  n_3D_timestep_start = n_3D_start;
end% if
if( n_3D_end == -1)
  N_3D_timestep = floor(TEND/DT);
  disp(['Looping up to rpinput-specified ' num2str(floor(TEND/DT))]);
else
  N_3D_timestep = n_3D_end;
  disp(['Looping from ' num2str(n_3D_start) ' to timestep ' num2str(n_3D_end) ' out of rpinput-specified ' num2str(floor(TEND/DT)) ' timesteps.']);
end% if
for n_3D_timestep = n_3D_timestep_start:DFQEBSLICE:N_3D_timestep,
% no info of fields etc for time step 0 (bean info only), thus put
% 0 at first
if( n_3D_timestep == 0)
  n_3D_timestep_str = sprintf('%.4d', DFQEBSLICE');
else
  n_3D_timestep_str = sprintf('%.4d', n_3D_timestep');
end% if
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
% establish quickpic dump plane (XZ or YZ)
qp_dump_plane = my_get_quickpic_dump_plane(datadir);
if (length (qp_version_suffix) == 0)
  disp('EA: unknown QuickPIC dump_plane');
  stop;
end% if
if( qp_dump_plane == 'XZ' )
  plane_ylabel = 'x [um]';
elseif( qp_dump_plane == 'YZ' )
  plane_ylabel = 'y [um]';
else
  plane_ylabel = '';
end% if
if(do_FALL)
  myfile = [datadir 'FEX-' qp_dump_plane '/FEX-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).FEX = double(my_read_hdf(myfile));
  myfile = [datadir 'FEY-' qp_dump_plane '/FEY-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).FEY = double(my_read_hdf(myfile));
  myfile = [datadir 'FEZ-' qp_dump_plane '/FEZ-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).FEZ = double(my_read_hdf(myfile));
  myfile = [datadir 'FBX-' qp_dump_plane '/FBX-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).FBX = double(my_read_hdf(myfile));
  myfile = [datadir 'FBY-' qp_dump_plane '/FBY-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).FBY = double(my_read_hdf(myfile));
  myfile = [datadir 'FBZ-' qp_dump_plane '/FBZ-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).FBZ = double(my_read_hdf(myfile));
end% if
if(do_QEB)
  myfile = [datadir 'QEB-' qp_dump_plane '/QEB-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).QEB = double(my_read_hdf(myfile));
end% if
if(do_QEB_3D)
  myfile = [datadir 'QEB/QEB_' n_3D_timestep_str qp_version_suffix];
  qp(n_3D_counter).QEB_3D = double(my_read_hdf(myfile));
end% if
if(do_QEP)
  if( strfind(qp_version_suffix, '.h5') )
    myfile = [datadir 'QEP1-' qp_dump_plane '/QEP1-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  else
    myfile = [datadir 'QEP01-' qp_dump_plane '/QEP01-' qp_dump_plane '_' n_3D_timestep_str qp_version_suffix];
  end% if
  qp(n_3D_counter).QEP = double(my_read_hdf(myfile));
end% if

%
% scale data to physical units
%
scale_x = Box_X / 2^INDX;
scale_y = Box_Y / 2^INDY;
scale_z = Box_Z / 2^INDZ;
offset_x0 = Box_X/2; % middle of box
offset_y0 = Box_Y/2; % middle of box
offset_z0 = Box_Z/2; % middle of box
scale_E = 1;
scale_B = 1;
scale_rho = 1;
% norm to cgs
omega_p = sqrt(n0*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
scale_E = scale_E * SI_em*1e3*SI_c*1e2*omega_p / (SI_e / 3.336e-10);
scale_B = scale_B * SI_em*1e3*SI_c*1e2*omega_p / (SI_e / 3.336e-10);
scale_rho = scale_rho * (SI_e / 3.336e-10) * n0*1e6;
% cgs to mks
scale_E = scale_E / ( 1e4 / SI_c);
scale_B = scale_B / 1e4;
scale_rho = scale_rho / (1e-5*SI_c);
% V/m to GV/m
scale_E = scale_E / 1e9;
% T to T [leave as is]
scale_B = scale_B;
% C/m^3 to N / cm^3
scale_rho = scale_rho/SI_e / 1e6;
% N / cm^3 to units of [np_0]
scale_rho = scale_rho / n0/1e6;
% update data
if(do_FALL)
  qp(n_3D_counter).FEX = qp(n_3D_counter).FEX * scale_E;
  qp(n_3D_counter).FEY = qp(n_3D_counter).FEY * scale_E;
  qp(n_3D_counter).FEZ = qp(n_3D_counter).FEZ * scale_E;
  qp(n_3D_counter).FBX = qp(n_3D_counter).FBX * scale_B;
  qp(n_3D_counter).FBY = qp(n_3D_counter).FBY * scale_B;
  qp(n_3D_counter).FBZ = qp(n_3D_counter).FBZ * scale_B;
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
n_3D_timestep_str_beam = sprintf('%.4d', n_3D_timestep');

% load phase space beam (merge all beam parts)
qp(n_3D_counter_beam).PP(n_beam).BEAM = my_get_quickpic_phasespace(datadir, n_beam_str, n_3D_timestep_str_beam); % auto-extracts beam part filenames exists
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
if( do_FEX )
  qp(n_3D_counter).FEX_max = max(max(qp(n_3D_counter).FEX(1:round(end*3/3), 1:end)));;
  qp(n_3D_counter).FEX_min = min(min(qp(n_3D_counter).FEX(1:round(end*3/3), 1:end)));
end% if
if( do_FBY )
  qp(n_3D_counter).FBY_max = max(max(qp(n_3D_counter).FBY(1:round(end*3/3), 1:end)));;
  qp(n_3D_counter).FBY_min = min(min(qp(n_3D_counter).FBY(1:round(end*3/3), 1:end)));
end% if
if( do_FALL )
  qp(n_3D_counter).FEX_max = max(max(qp(n_3D_counter).FEX(1:round(end*3/3), 1:end)));;
  qp(n_3D_counter).FEX_min = min(min(qp(n_3D_counter).FEX(1:round(end*3/3), 1:end)));
  qp(n_3D_counter).FEY_max = max(max(qp(n_3D_counter).FEY(1:round(end*3/3), 1:end)));;
  qp(n_3D_counter).FEY_min = min(min(qp(n_3D_counter).FEY(1:round(end*3/3), 1:end)));
  qp(n_3D_counter).FEZ_max = max(max(qp(n_3D_counter).FEZ(1:round(end*3/3), 1:end)));;
  qp(n_3D_counter).FEZ_min = min(min(qp(n_3D_counter).FEZ(1:round(end*3/3), 1:end)));
  qp(n_3D_counter).FBX_max = max(max(qp(n_3D_counter).FBX(1:round(end*3/3), 1:end)));;
  qp(n_3D_counter).FBX_min = min(min(qp(n_3D_counter).FBX(1:round(end*3/3), 1:end)));
  qp(n_3D_counter).FBY_max = max(max(qp(n_3D_counter).FBY(1:round(end*3/3), 1:end)));;
  qp(n_3D_counter).FBY_min = min(min(qp(n_3D_counter).FBY(1:round(end*3/3), 1:end)));
  qp(n_3D_counter).FBZ_max = max(max(qp(n_3D_counter).FBZ(1:round(end*3/3), 1:end)));;
  qp(n_3D_counter).FBZ_min = min(min(qp(n_3D_counter).FBZ(1:round(end*3/3), 1:end)));
end% if
qp(n_3D_counter).QEB_max =  max(max(qp(n_3D_counter).QEB));
qp(n_3D_counter).PP(n_beam).mean_E = mean(qp(n_3D_counter_beam).PP(n_beam).BEAM(:,6));
% store slice transverse data
if(n_3D_counter == 1)
  [slice.mean_x, slice.sigma_x, slice.z(n_beam, :), slice.N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 1);
else
  % ensure same slicing as beam evolves through plasma, in order to compare apples to apples
  pp = qp(n_3D_counter_beam).PP(n_beam).BEAM;
  [Slice.mean_x, slice.sigma_x, slice.z(n_beam, :), slice.N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 1, slice.z(n_beam, :));
end% if
[slice.mean_y, slice.sigma_y, slice.z(n_beam, :), slice.N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 2);
qp(n_3D_counter).PP(n_beam).slice.z = slice.z(n_beam, :);
qp(n_3D_counter).PP(n_beam).slice.N_z = slice.N_z;
qp(n_3D_counter).PP(n_beam).slice.mean_x = slice.mean_x;
qp(n_3D_counter).PP(n_beam).slice.sigma_x = slice.sigma_x;
qp(n_3D_counter).PP(n_beam).slice.mean_y = slice.mean_y;
qp(n_3D_counter).PP(n_beam).slice.sigma_y = slice.sigma_y;
% store slice angle data
[slice.mean_xp, slice.sigma_xp, slice.z(n_beam, :), slice.N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 4, slice.z(n_beam, :));
[slice.mean_yp, slice.sigma_yp, slice.z(n_beam, :), slice.N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 5, slice.z(n_beam, :));
qp(n_3D_counter).PP(n_beam).slice.mean_xp = slice.mean_xp;
qp(n_3D_counter).PP(n_beam).slice.sigma_xp = slice.sigma_xp;
qp(n_3D_counter).PP(n_beam).slice.mean_yp = slice.mean_yp;
qp(n_3D_counter).PP(n_beam).slice.sigma_yp = slice.sigma_yp;
% store slice energy data
if(n_3D_counter == 1)
  [slice.mean, slice.sigma, slice.z(n_beam, :), slice.N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 6);
else
  [slice.mean, slice.sigma, slice.z(n_beam, :), slice.N_z] = my_get_slice_var(qp(n_3D_counter_beam).PP(n_beam).BEAM, 3, 6, slice.z(n_beam, :));
end% if
qp(n_3D_counter).PP(n_beam).slice.z = slice.z(n_beam, :);
qp(n_3D_counter).PP(n_beam).slice.N_z = slice.N_z;
qp(n_3D_counter).PP(n_beam).slice.mean_E = slice.mean;
qp(n_3D_counter).PP(n_beam).slice.sigma_E = slice.sigma;
% store slice twiss data
[emnx,emny,betax,alphax,betay,alphay] = my_get_slice_twiss(qp(n_3D_counter_beam).PP(n_beam).BEAM, slice.z(n_beam, :));
qp(n_3D_counter).PP(n_beam).slice.emnx = emnx;
qp(n_3D_counter).PP(n_beam).slice.emny = emny;
qp(n_3D_counter).PP(n_beam).slice.betax = betax;
qp(n_3D_counter).PP(n_beam).slice.alphax = alphax;
qp(n_3D_counter).PP(n_beam).slice.betay = betay;
qp(n_3D_counter).PP(n_beam).slice.alphay = alphay;
% store slice gauss fit data (time consuming)
if( do_ana_slice__gauss_fit )
  [gaussx, gaussxp, gaussy, gaussyp] = my_get_slice_gauss(qp(n_3D_counter_beam).PP(n_beam).BEAM, slice.z(n_beam, :));
  qp(n_3D_counter).PP(n_beam).slice.gaussx = gaussx;
  qp(n_3D_counter).PP(n_beam).slice.gaussxp = gaussxp;
  qp(n_3D_counter).PP(n_beam).slice.gaussy = gaussy;
  qp(n_3D_counter).PP(n_beam).slice.gaussyp = gaussyp;
  mysgx(n_3D_counter) = gaussx(30);
  mysgy(n_3D_counter) = gaussy(30);
end% if

mybx(n_3D_counter) = betax(30);
myby(n_3D_counter) = betay(30);
myex(n_3D_counter) = emnx(30);
myey(n_3D_counter) = emny(30);
mybx(n_3D_counter) = betax(30);
myby(n_3D_counter) = betay(30);
myssx(n_3D_counter) = slice.sigma_x(30);
myssy(n_3D_counter) = slice.sigma_y(30);
if( do_FULL_beam_phase_space_plot && do_plot )
  pp = qp(1).PP(n_beam).BEAM;
  [Y, I] = sort(pp(:,3));
  pp = pp(I, :);
  my_ana_beam(pp(end*0/2+1:end,:), [1 1 1 0 1 0 0]);
  pause;
end; % if
if( do_SLICE_beam_phase_space_plot && do_plot )
  pp = qp(1).PP(n_beam).BEAM;
  [Y, I] = sort(pp(:,3));
  pp = pp(I, :);
  n_min = max(find(Y < qp(n_3D_counter).PP(n_beam).slice.z(n_slice_plot)));
  n_max = max(find(Y < qp(n_3D_counter).PP(n_beam).slice.z(n_slice_plot+1)));
  pp_slice = pp(n_min:n_max, :);
  my_ana_beam(pp_slice, [1 1 1 0 1 0 0]);
  pause;
end; % if
% intital charge and gamma, per beam
qp(n_3D_counter).PP(n_beam).N0 = qp_N(n_beam);
qp(n_3D_counter).PP(n_beam).gamma0 = qp_Gamma(n_beam);
qp(n_3D_counter).E_dec_max =  max(max(qp(n_3D_counter).FEZ(1:round(end*1/2), 1:end))); % peak dec field
% fraction of accelerated and decelerated charge
[Y, I] = sort(qp(1).PP(n_beam).BEAM(:,6));
pp_E = qp(1).PP(n_beam).BEAM(I, 6);
E0 = qp(n_3D_counter).PP(n_beam).gamma0*.511e-3;
sigma_E_E0 = 0.04; %  
E_dec = E0*(1-sigma_E_E0);
E_acc = E0*(1+sigma_E_E0);
s_timestep = (SI_c/omega_p) * DT * n_3D_timestep; % propagation length for this timestep
E_acc_T2 = E0 + 2*qp(n_3D_counter).E_dec_max*s_timestep; % Energy for particles experiencing a T>=2
n_dec_max = max(find(pp_E<E_dec));
n_acc_min = min(find(pp_E>E_acc));
n_acc_T2_min = min(find(pp_E>E_acc_T2));
if( length(n_acc_min) > 0)
  n_acc_frac = 1-n_acc_min / length(pp_E);
else
  n_acc_frac = 0;
end% if
if( length(n_acc_T2_min) > 0)
  n_acc_T2_frac = 1-n_acc_T2_min / length(pp_E);
else
  n_acc_T2_frac = 0;
end% if
if( length(n_dec_max) > 0)
  n_dec_frac = n_dec_max / length(pp_E);
else
  n_dec_frac = 0;
end% if
qp(n_3D_counter).PP(n_beam).n_acc_frac = n_acc_frac;
qp(n_3D_counter).PP(n_beam).n_acc_T2_frac = n_acc_T2_frac;
qp(n_3D_counter).PP(n_beam).n_dec_frac = n_dec_frac;

end% for each beam

% generic values
qp(n_3D_counter).timestamp = now(); % timestamp data file
qp(n_3D_counter).n0 = n0;
qp(n_3D_counter).neutral_gas = neutral_gas; 
s_timestep = (SI_c/omega_p) * DT * n_3D_timestep; % propagation length for this timestep
qp(n_3D_counter).s_timestep = s_timestep;
disp(' ' ); 
display(['Prop. length s=' num2str(s_timestep*100, 3) ' [cm]. Time steps in units of DT: ' n_3D_timestep_str]);

%
% calculating hosing and ctm quantities
%
if(do_calc_hose_ctm)
  do_beam_1_only = 1;
  if( do_beam_1_only == 1)
    pp = qp(n_3D_counter_beam).PP(1).BEAM;
  else
    pp = qp(n_3D_counter_beam).PP(n_beam).BEAM;
  end% if
[Y, I] = sort(pp(:,3));
pp = pp(I, :);
%zn_ion = min(find(sum(qp(n_3D_counter).QEP', 1) > 0));
% Not working for highly sloshing beam
%[val, xn_ion] = max(max(qp(n_3D_counter).FEZ(1:round(end*2/3), :)));
%zn_ion = min(find(qp(n_3D_counter).FEZ(1:round(end*2/3), xn_ion) > max(qp(n_3D_counter).FEZ(1:round(end*2/3), xn_ion))*1/2 ));
% Not working for highly sloshing beam
dens_peak = max(max(qp(n_3D_counter).QEB));
zn_ion = min(find(max(qp(n_3D_counter).QEB') >  1/2*dens_peak));
dz0_offset = 0; % [um], start 
dz0_length = 10; % [um], start 
z_ion_head = zn_ion*scale_z-offset_z0 + dz0_offset;
% if first timestep (no plasma interaction), set inital z_ion_head in
% determinstic way
if( n_3D_timestep == 0)
  sgh0=8.0;
  z_ion_head = pp( round(size(pp,1)*sgh0/100),3 ) + dz0_offset;
end% if
z1_ion_head = zn_ion*scale_z-offset_z0 + (dz0_offset+dz0_length);
nn0 = max(find( pp(:,3) < z_ion_head));
nn1 = max(find( pp(:,3) < z1_ion_head));
z_ion_head;
x_ion_head =  mean(pp(nn0:nn1,1));
xp_ion_head =  mean(pp(nn0:nn1,4));
y_ion_head =  mean(pp(nn0:nn1,2));
yp_ion_head =  mean(pp(nn0:nn1,5));
% if first timestep (no plasma interaction), set initalion_head in
% determinstic way
if( n_3D_timestep == 0)
  x_ion_head =  0;
  y_ion_head =  0;
  xp_ion_head =  0;
  yp_ion_head =  0;
end% if

% fix tail particles to be tracked at first time step
if( n_3D_counter == 1)
  % at this percentage of the beam, and keep same n_part throughout?
  sgh0 = 84.2;
  sgh1 = 86.3;
  n_tail_start_0 = round(size(pp,1)*sgh0/100);
  n_tail_end_0 = round(size(pp,1)*sgh1/100);
  % since head erodes -> lost particles, better to keep fixed z
  z_tail_start = pp(n_tail_start_0,3);
  z_tail_end = pp(n_tail_end_0,3);
  % or simply, set fixed z directly : at n initial sigmas after mean of beam
  n_sigma = 1;
  qp(1).PP(n_beam).hose_calc.n_sigma = n_sigma;
  z_tail_start = qp(1).PP(1).mean_z + (n_sigma-0.1)*qp(1).PP(1).sigma_z;
  z_tail_end = qp(1).PP(1).mean_z + (n_sigma+0.1)*qp(1).PP(1).sigma_z;
end% if

n_tail_start = max(find(pp(:,3) < z_tail_start));
n_tail_end = max(find(pp(:,3) < z_tail_end));
%n_tail_start = n_tail_start_0;
%n_tail_end = n_tail_end_0;
x_tail =  mean(pp(n_tail_start:n_tail_end,1));
y_tail =  mean(pp(n_tail_start:n_tail_end,2));
qp(n_3D_counter).PP(n_beam).x_tail = x_tail;
qp(n_3D_counter).PP(n_beam).y_tail = y_tail;
qp(n_3D_counter).PP(n_beam).x_head = mean(pp(1:round(size(pp,1)*1/100),1));
qp(n_3D_counter).PP(n_beam).y_head = mean(pp(1:round(size(pp,2)*1/100),2));
qp(n_3D_counter).PP(n_beam).x_ion_head = x_ion_head;
qp(n_3D_counter).PP(n_beam).y_ion_head = y_ion_head;
qp(n_3D_counter).PP(n_beam).xp_ion_head = xp_ion_head;
qp(n_3D_counter).PP(n_beam).yp_ion_head = yp_ion_head;
qp(n_3D_counter).PP(n_beam).z_ion_head = z_ion_head;
 
nn0h = round(size(pp,1)*sgh0/100);
nn1h = round(size(pp,1)*sgh1/100);
do_thisplot = 0;
if( do_thisplot )
  plot(pp(:,3), pp(:,1), '.');
  hold on;
  plot(pp(nn0:nn1,3), pp(nn0:nn1,1), '.r');
  plot(pp(nn0h:nn1h,3), pp(nn0h:nn1h,1), '.g');
  hold off;
  xlabel('z [um]');
  ylabel('x [um]');
  grid on;
  pause
end% if

%
% calc' bubble and focusing field calcs
%
  ZZ = (1:size(qp(n_3D_counter).QEB,1)) * scale_z - offset_z0;
  XX = (1:size(qp(n_3D_counter).QEB,2)) * scale_x - offset_x0;
  YY = (1:size(qp(n_3D_counter).QEB,2)) * scale_y - offset_y0;
  if( qp_dump_plane == 'XZ' )
    RR = XX;
    EE = qp(n_3D_counter).FEX; % [GV/m]
    BB = -SI_c*qp(n_3D_counter).FBY/1e9; % [GV/m]
  elseif( qp_dump_plane == 'YZ' )
    RR = YY;
    EE = qp(n_3D_counter).FEY; % [GV/m]
    BB = +SI_c*qp(n_3D_counter).FBX/1e9; % [GV/m]
  end% if
  % calc EB slope either at fixed number, or at erosion front
  if( strcmp(translopeEB_dz, 'ion_head'))
    % at erosion front
    translopeEB_dz_forcalc = (z1_ion_head + z_ion_head)/2;
    x_beam = x_ion_head;
  elseif( strcmp(translopeEB_dz, 'hose_tail'))
    % at hose calc point
    translopeEB_dz_forcalc = (z_tail_end + z_tail_start)/2;
    x_beam = x_tail;
  else
    % fixed number
    translopeEB_dz_forcalc = translopeEB_dz;
    x_beam = 0; 
  end% if
  d_INDZ = round(translopeEB_dz_forcalc/scale_z);
  EE_BB =  EE((2^INDZ/2)+1+d_INDZ, :) + BB((2^INDZ/2)+1+d_INDZ, :);
  % robust bubble radius: deduce bubble radius from field itself (assume good limits)
  [VL, IXp] = max(EE_BB);
  [VL, IXm] = min(EE_BB);
  R0 = (IXm + my_find_zero_cross(EE_BB(IXm:IXp), 1) - 1) * scale_x;
  Rp = IXp*scale_x;
  Rm = IXm*scale_x;
  R_bubble_calc = (Rp-Rm)/2;
  f_linear = 0.50;
  if( isstr( translopeEB_dz ))
    % at erosion front
    translopeEB_Rp = +f_linear*abs(R_bubble_calc) + x_beam;
    translopeEB_Rm = -f_linear*abs(R_bubble_calc) + x_beam;
  else
    % fixed number
    translopeEB_Rp = +f_linear*abs(R_bubble_calc); % assume linear with 60% from max
    translopeEB_Rm = -f_linear*abs(R_bubble_calc); % assume linear with 60% from max
  end% if
  x0 = (length(EE( (2^INDZ/2)+1+d_INDZ))/2 + 1) * scale_x + x_beam;
  xp = (length(EE( (2^INDZ/2)+1+d_INDZ))/2 + 1) * scale_x + translopeEB_Rp;
  xm = (length(EE( (2^INDZ/2)+1+d_INDZ))/2 + 1) * scale_x + translopeEB_Rm;
  E_R0 = EE( (2^INDZ/2)+1+d_INDZ, end/2+1);
  E_Rp = EE( (2^INDZ/2)+1+d_INDZ, end/2+1 + round(translopeEB_Rp/scale_x) );
  E_Rm = EE( (2^INDZ/2)+1+d_INDZ, end/2+1 + round(translopeEB_Rm/scale_x) );
  B_R0 = BB( (2^INDZ/2)+1+d_INDZ, end/2+1);
  B_Rp = BB( (2^INDZ/2)+1+d_INDZ, end/2+1 + round(translopeEB_Rp/scale_x) );
  B_Rm = BB( (2^INDZ/2)+1+d_INDZ, end/2+1 + round(translopeEB_Rm/scale_x) );
  P = polyfit( [xm x0 xp], [E_Rm+B_Rm E_R0+B_R0 E_Rp+B_Rp], 1);
  % if we look at field at beam motion, let's mark where the beam <x> is
  %    by putting the center of the fit line here.
  if( isstr( translopeEB_dz ))
    dx = x_beam - x0;
    dy = 0 - (P(2) + P(1)*x0);
  else
    dx = 0;
    dy = 0;
  end% if
  % for which x is the field zero?
  [dummy, indx_min] = min(abs(EE_BB((end/2+1 + round(translopeEB_Rm/scale_x)) : (end/2+1 + round(translopeEB_Rp/scale_x)))));
  indx_min = indx_min + length(EE_BB)/2 + round(translopeEB_Rm/scale_x);
  x_zero_EcB = RR(indx_min);
  % we are not fully at zero due to granularit of grid, so store
  EcB_zero =  EE_BB(indx_min);
  qp(n_3D_counter).PP(n_beam).x_zero_EcB = x_zero_EcB;
  qp(n_3D_counter).PP(n_beam).EcB_zero = EcB_zero;
  % what is the field at the beam x?
  if( strcmp(translopeEB_dz, 'ion_head'))
    % at erosion front
    dx_zeros_beam_field = x_beam - x_zero_EcB;
  elseif( strcmp(translopeEB_dz, 'hose_tail'))
    % at erosion front
    dx_zeros_beam_field = x_tail - x_zero_EcB;
  else
    % fixed number
    dx_zeros_beam_field = 0;
  end% if
  EcB_at_beam = EE_BB(round(indx_min + dx_zeros_beam_field/scale_x));
  qp(n_3D_counter).PP(n_beam).EcB_at_beam = EcB_at_beam;
  qp(n_3D_counter).PP(n_beam).x_beam_EcB = x_beam;
  % what is the field at the hosing calc point? [run again to see]
  
  % theoretical value, to compare (from Esaray, I. Blumenfeld)
  F_over_e_m_theory = (1/2)*(SI_em)*SI_c^2*k_p^2 / SI_e / SI_c; % F/e/L in [T/m] (mult by SI_c/1e9 for [T/m])
  dB_dx = P(1)/SI_c * 1e15; % [T/M]


end% if

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
I_peak = qp_N(1)*SI_e *SI_c/sqrt(2*pi*(qp(n_3D_counter).PP(n_beam).sigma_z/1e6)^2);
n_dens_bunch = qp_N(1) / ( (2*pi)^(3/2)*qp(n_3D_counter).PP(n_beam).sigma_x/1e6*qp(n_3D_counter).PP(n_beam).sigma_y/1e6*qp(n_3D_counter).PP(n_beam).sigma_z/1e6 ) / 1e6; % cm^-3
lambda_betatron = sqrt(2*qp_Gamma(1))*(2*pi) / k_p;
Lambda_nonlin = 2.5*(qp_N(1)/2e10)*(20/qp(n_3D_counter).PP(n_beam).sigma_x);
crit_blow_out = n_dens_bunch / (n0 * (1+4/(k_p*qp(n_3D_counter).PP(n_beam).sigma_z/1e6)^2) );
V_erosion_PI = sqrt( qp(n_3D_counter).PP(n_beam).emitt(1)*k_p/(qp(n_3D_counter).PP(n_beam).mean_E/.511e-3)^(3/2)/(I_peak/17e3) );
E0_wb = SI_em*SI_c*omega_p / SI_e / 1e9; % Wave-breaking field, GV/m
Ez_lin_IB = 4*pi*E0_wb * I_peak / 17e3;
Ez_lin = -100*(qp_N(1)/2e10)*(20/qp(n_3D_counter).PP(n_beam).sigma_z)^2*log(sqrt(2.5e17*1e6 / n0/1e6 * 10 / qp(n_3D_counter).PP(n_beam).sigma_x)); % lin. regime (NJP), GV/m
Ez_sel = 40*(qp_N(1)/2e10)*(15/qp(n_3D_counter).PP(n_beam).sigma_x)*(20/qp(n_3D_counter).PP(n_beam).sigma_z)^2;
Ez_hog = 0.244*(qp_N(1)/2e10)*(600/qp(n_3D_counter).PP(n_beam).sigma_z)^2;
% hose formulae
omega_p = k_p / sqrt(2);
k_beta = k_p / sqrt(2*qp_Gamma(1));
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
  subplot(2,3,2);

% TEMP IPAC'12
%figure1 = figure('Color',[1 1 1]);
%subplot(1,1,1);
%set(0,'defaultaxesfontsize',24);
%qp(n_3D_counter).QEB = flipud(qp(n_3D_counter).QEB);
%qp(n_3D_counter).QEP = flipud(qp(n_3D_counter).QEP);
  
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

  map_max = size(eval(mycolormap),1);
%  mycolormap2 = [gray(map_max); jet(map_max)];
%  colormap(mycolormap2);
%  pcolor(ZZ,XX,qp(n_3D_counter).QEB');
  pcolor(ZZ,XX,((qp(n_3D_counter).QEB'*((map_max-1)/max(max(qp(n_3D_counter).QEB))) + qp(n_3D_counter).QEP'*((map_max/16-1)/mean(mean(qp(n_3D_counter).QEP))))));  % with plasma
  caxis([0 (map_max)]);
%  pcolor(ZZ,XX,qp(n_3D_counter).QEP');  % ONLY plasma
%  pcolor(ZZ,XX,qp(n_3D_counter).QEB');  % only beam
  h_c = colorbar('NorthOutside');
  %caxis([0.0 0.1]);
  if(do_fixed_axes)
    caxis([0 100]); % NJP
    caxis([0 30]);
    caxis([0 5]); % Dx-3
    caxis([0 5]); % Dx-0 mat
    caxis([-1.0 1.1]); % positrons
  end% if
%    caxis([0 10]);
%    caxis([0 4]); % positrons

  % IPAC'12
%  h_c = colorbar('EastOutside');
%    caxis([0 5]); % Dx-0 mat
%%    axis([-100 50 -100 100])
%    axis([-50 100 -100 100])
%title(['s=' num2str(s_timestep*100, 3) ' [cm].   Time step: ' num2str(n_3D_timestep) ' [DT]']);

shading('flat');
  xlabel('z [um]')
  ylabel( plane_ylabel)
  %set(get(h_c,'ylabel'),'String', 'n_b  [n_p]');
  %set(get(h_c,'ylabel'),'String', 'n_b  [C/cm^3]');
%  set(get(h_c,'ylabel'),'String', 'n_b  [n_p_0]');
  set(get(h_c,'ylabel'),'String', 'n_e  [n_p_0]');  % with plasma

% IPAC'12  
%  stop
  
  % add line plot on axis
  subplot(2,3,5);
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
    axis([min(ZZ)  max(ZZ)  0 20]);
    %axis([min(ZZ)  max(ZZ)  0 QEB_max_0]); % Dx03
    axis([min(ZZ)  max(ZZ)  -10 10]); % positrons
  end% if
  QEB_max =  max(max(qp(n_3D_counter).QEB));
  title(['n_{b,max} = ' num2str(QEB_max, '%.2f') ' n_p_0']);
%  pause;
end% if plot


%
%  plot QEP
%
if(do_QEP && 0)
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
  ylabel( plane_ylabel)
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
  subplot(2,3,3);


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
  ylabel( plane_ylabel)
  %set(get(h_c,'ylabel'),'String', 'n_b  [n_p]');
  %set(get(h_c,'ylabel'),'String', 'n_b  [C/cm^3]');
  set(get(h_c,'ylabel'),'String', 'E_z [GV/m]');

  % add line plot on axis
  subplot(2,3,6);

      % TEMP IPAC'12
%qp(n_3D_counter).FEZ = flipud(qp(n_3D_counter).FEZ);
%figure1 = figure('Color',[1 1 1]);
%set(0,'defaultaxesfontsize',24);
%
  hh = plot(ZZ, qp(n_3D_counter).FEZ(:,(2^INDX/2)+1));
  grid on;
  xlabel('z [um]')
  ylabel('E_z [GV/m]')
  qp_FEZ_max = max(qp_FEZ_max, max(qp(n_3D_counter).FEZ(:,(2^INDX/2)+1)));
  qp_FEZ_min = min(qp_FEZ_min, min(qp(n_3D_counter).FEZ(:,(2^INDX/2)+1)));
  myaxis = axis;
  axis([min(ZZ)  max(ZZ)  qp_FEZ_min  qp_FEZ_max+eps]);
  FEZ_max_0 = max(max(qp(1).FEZ(1:round(end*1/3), 1:end)));
  if(do_fixed_axes)
    axis([min(ZZ)  max(ZZ) -70 70]); % NJP
%    axis([min(ZZ)  max(ZZ) -50 50]);
%    axis([min(ZZ)  max(ZZ) -FEZ_max_0 FEZ_max_0]); % Dx-3
  end% if
  FEZ_max = max(max(qp(n_3D_counter).FEZ(1:round(end*1/2), 1:end)));
  FEZ_min = min(min(qp(n_3D_counter).FEZ(1:round(end*1/2), 1:end)));
  title(['E_{z,max} = ' num2str(FEZ_max, '%.1f') 'GV/m']);
  %title(['E_{z,max} = ' num2str(FEZ_max, '%.1f') 'GV/m'  ', ' 'E_{z,min} = ' num2str(FEZ_min, '%.1f') 'GV/m']);

  % TEMP IPAC'12
%  set(hh,'LineWidth',5)
%    axis([-50 100 -120 120])
%  pause;

end% if plot


%
% TEMP QEB + FEZ
%
if(do_QEB_FEZ)
  n_zc = my_find_zero_cross(qp(n_3D_counter).FEZ(:,(2^INDX/2)+1), -1);
  n_frac_acc = sum(qp(n_3D_counter).QEB(n_zc:end,:)) / sum(qp(n_3D_counter).QEB(:,:));
  
  clf;
  hh = plot(ZZ, qp(n_3D_counter).QEB(:,(2^INDX/2)+1)*10, '-b'); % middle line
  set(hh, 'linewidth', 4)
  hold on;
  hh = plot(ZZ, qp(n_3D_counter).FEZ(:,(2^INDX/2)+1), '-r');
  set(hh, 'linewidth', 4)
hold off;
xlabel('z [um]');
ylabel('10*n_{eb}/n_{0} [-],   E_{z} [GV/m]');
legend('beam charge dens [n_{0}/10]', 'plasma wake field [GV/m]');
title(['\sigma_{z} [um]=20,' 'n_{0} [/1e17xcm^3]=' num2str(qp(n_3D_counter).n0/1e17) '  : frac Q in ACC field [%] :' num2str(n_frac_acc*100,2)]);
grid on;
axis([-200 200 -100 100]);
  stop
end% if



%
%  plot FEX
%
if(do_FEX)
  subplot(2,3,3);
  colormap(mycolormap);
  ZZ = (1:size(qp(n_3D_counter).FEX,1)) * scale_z - offset_z0;
  XX = (1:size(qp(n_3D_counter).FEX,2)) * scale_x - offset_x0;
  YY = (1:size(qp(n_3D_counter).FEX,2)) * scale_y - offset_y0;
  c_min_new = (min( min(min(qp(n_3D_counter).FEX)), min(min(qp(n_3D_counter).FEX)) ));
  c_max_new = (max( max(max(qp(n_3D_counter).FEX)), max(max(qp(n_3D_counter).FEX)) ));
  % don't jump up and down as time evolve
  if(c_min_new < c_min)
    c_min = c_min_new;
  end% if
  if(c_max_new > c_max)
    c_max = c_max_new;
  end% if

  pcolor(ZZ,XX,qp(n_3D_counter).FEX');
  h_c = colorbar('NorthOutside');
%  caxis([c_min c_max]);
  if(do_fixed_axes)
    caxis([-70 70]); % NJP
    caxis([-30 30]);
%    caxis([-5 5]); % Dx-3
%    caxis([-5 5]); % positrons
  end% if
  shading('flat');
  xlabel('z [um]')
  ylabel( plane_ylabel)
  %set(get(h_c,'ylabel'),'String', 'n_b  [n_p]');
  %set(get(h_c,'ylabel'),'String', 'n_b  [C/cm^3]');
  set(get(h_c,'ylabel'),'String', 'E_x [GV/m]');

  % add line plot on axis
  subplot(2,3,6);
  plot(ZZ, qp(n_3D_counter).FEX(:,(2^INDX/2)+1));
  grid on;
  xlabel('z [um]')
  ylabel('E_x [GV/m]')
  qp_FEX_max = max(qp_FEX_max, max(qp(n_3D_counter).FEX(:,(2^INDX/2)+1)));
  qp_FEX_min = min(qp_FEX_min, min(qp(n_3D_counter).FEX(:,(2^INDX/2)+1)));
  myaxis = axis;
  axis([min(ZZ)  max(ZZ)  qp_FEX_min  qp_FEX_max+eps]);
  FEX_max_0 = max(max(qp(1).FEX(1:round(end*2/3), 1:end)));
  if(do_fixed_axes)
    axis([min(ZZ)  max(ZZ) -70 70]); % NJP
%    axis([min(ZZ)  max(ZZ) -50 50]);
%    axis([min(ZZ)  max(ZZ) -FEX_max_0 FEX_max_0]); % Dx-3
  end% if
  FEX_max = max(max(qp(n_3D_counter).FEX(1:round(end*2/3), 1:end)));
  FEX_min = min(min(qp(n_3D_counter).FEX(1:round(end*2/3), 1:end)));
  title(['E_{x,max} = ' num2str(FEX_max, '%.1f') 'GV/m']);
  %title(['E_{x,max} = ' num2str(FEX_max, '%.1f') 'GV/m'  ', ' 'E_{x,min} = ' num2str(FEX_min, '%.1f') 'GV/m']);

%  pause;
end% if plot



%
%  plot FEY
%
if(do_FEY)
  subplot(2,3,3);
  colormap(mycolormap);
  ZZ = (1:size(qp(n_3D_counter).FEY,1)) * scale_z - offset_z0;
  XX = (1:size(qp(n_3D_counter).FEY,2)) * scale_x - offset_x0;
  YY = (1:size(qp(n_3D_counter).FEY,2)) * scale_y - offset_y0;
  c_min_new = (min( min(min(qp(n_3D_counter).FEY)), min(min(qp(n_3D_counter).FEY)) ));
  c_max_new = (max( max(max(qp(n_3D_counter).FEY)), max(max(qp(n_3D_counter).FEY)) ));
  % don't jump up and down as time evolve
  if(c_min_new < c_min)
    c_min = c_min_new;
  end% if
  if(c_max_new > c_max)
    c_max = c_max_new;
  end% if

  pcolor(ZZ,XX,qp(n_3D_counter).FEY');
  h_c = colorbar('NorthOutside');
%  caxis([c_min c_max]);
  if(do_fixed_axes)
    caxis([-70 70]); % NJP
    caxis([-30 30]);
%    caxis([-5 5]); % Dx-3
%    caxis([-5 5]); % positrons
  end% if
  shading('flat');
  xlabel('z [um]')
  ylabel( plane_ylabel)
  %set(get(h_c,'ylabel'),'String', 'n_b  [n_p]');
  %set(get(h_c,'ylabel'),'String', 'n_b  [C/cm^3]');
  set(get(h_c,'ylabel'),'String', 'E_x [GV/m]');

  % add line plot on axis
  subplot(2,3,6);
  plot(ZZ, qp(n_3D_counter).FEY(:,(2^INDX/2)+1));
  grid on;
  xlabel('z [um]')
  ylabel('E_x [GV/m]')
  qp_FEY_max = max(qp_FEY_max, max(qp(n_3D_counter).FEY(:,(2^INDX/2)+1)));
  qp_FEY_min = min(qp_FEY_min, min(qp(n_3D_counter).FEY(:,(2^INDX/2)+1)));
  myaxis = axis;
  axis([min(ZZ)  max(ZZ)  qp_FEY_min  qp_FEY_max+eps]);
  FEY_max_0 = max(max(qp(1).FEY(1:round(end*2/3), 1:end)));
  if(do_fixed_axes)
    axis([min(ZZ)  max(ZZ) -70 70]); % NJP
%    axis([min(ZZ)  max(ZZ) -50 50]);
%    axis([min(ZZ)  max(ZZ) -FEY_max_0 FEY_max_0]); % Dx-3
  end% if
  FEY_max = max(max(qp(n_3D_counter).FEY(1:round(end*2/3), 1:end)));
  FEY_min = min(min(qp(n_3D_counter).FEY(1:round(end*2/3), 1:end)));
  title(['E_{x,max} = ' num2str(FEY_max, '%.1f') 'GV/m']);
  %title(['E_{x,max} = ' num2str(FEY_max, '%.1f') 'GV/m'  ', ' 'E_{x,min} = ' num2str(FEY_min, '%.1f') 'GV/m']);

%  pause;
end% if plot





%
%  plot and estimate trans slope 
%
if(do_translopeEB)
  % add line plot on axis
  subplot(2,3,3);
  plot(RR, EE_BB, '-r');
  hold on;
  set(hh,'MarkerSize', 30)
  myaxis = axis;
    axis([x_beam-2*f_linear*abs(R_bubble_calc) x_beam+2*f_linear*abs(R_bubble_calc) myaxis(3) myaxis(4)]);
    grid on;
  xlabel([plane_ylabel ' @ z=' num2str(translopeEB_dz_forcalc, 2) ' um'])
  ylabel('E - cB [GV/m]')
  hh = plot([xm x0 xp]+dx, (P(2) + P(1)*([xm x0 xp]))+dy, '-xk')';
  set(hh,'LineWidth',3)
  % plot beam zero
  hh = plot(x0+dx, (P(2) + P(1)*([x0]))+dy, '+k')';
  set(hh,'MarkerSize', 30)
  % plot field zero
  hh = plot(RR(indx_min), EE_BB(indx_min), '+r');
  set(hh,'MarkerSize', 30);
end% if




%
%  plot E x B
%
if(do_EcrossB)
  % calc E x B
  disp('...calcing E x B; please wait ~ 10 sec');
  ExB = my_calc_cross(qp(n_3D_counter).FEX, qp(n_3D_counter).FEY, qp(n_3D_counter).FEZ, qp(n_3D_counter).FBX, qp(n_3D_counter).FBY, qp(n_3D_counter).FBZ);
  %qp(n_3D_counter).ExB = my_structmat2mat(ExB, 5); % |E x B|
  qp(n_3D_counter).v_ExB = my_structmat2mat(ExB, 7); % |v|
                            
  subplot(2,3,3);
%  subplot(2,1,1);
  colormap(mycolormap);
  ZZ = (1:size(qp(n_3D_counter).v_ExB,1)) * scale_z - offset_z0;
  XX = (1:size(qp(n_3D_counter).v_ExB,2)) * scale_x - offset_x0;
  YY = (1:size(qp(n_3D_counter).v_ExB,2)) * scale_y - offset_y0;
  c_min_new = (min( min(min(qp(n_3D_counter).v_ExB)), min(min(qp(n_3D_counter).v_ExB)) ));
  c_max_new = (max( max(max(qp(n_3D_counter).v_ExB)), max(max(qp(n_3D_counter).v_ExB)) ));
  % don't jump up and down as time evolve
  if(c_min_new < c_min)
    c_min = c_min_new;
  end% if
  if(c_max_new > c_max)
    c_max = c_max_new;
  end% if

  % cell # corresponding to beam centroid
  qp(n_3D_counter).PP(n_beam).centr_cell_x = max(find(XX < qp(n_3D_counter).PP(n_beam).mean_x));
  qp(n_3D_counter).PP(n_beam).centr_cell_y = max(find(YY < qp(n_3D_counter).PP(n_beam).mean_y));
  qp(n_3D_counter).PP(n_beam).centr_cell_z = max(find(ZZ < qp(n_3D_counter).PP(n_beam).mean_z));
  qp(n_3D_counter).PP(n_beam).centr_cell_mag_v = ExB( qp(n_3D_counter).PP(n_beam).centr_cell_x, qp(n_3D_counter).PP(n_beam).centr_cell_z).mag_v ;
  
  pcolor(ZZ,XX,qp(n_3D_counter).v_ExB');
  h_c = colorbar('NorthOutside');
%  caxis([c_min c_max]);
  if(do_fixed_axes)
    caxis([0 50]); % NJP
    caxis([0 5]);
%    caxis([-5 5]); % Dx-3
%    caxis([-5 5]); % positrons
  end% if
    caxis([0 5]);
  shading('flat');
  Xlabel('z [um]')
  ylabel( plane_ylabel)
  %set(get(h_c,'ylabel'),'String', 'n_b  [n_p]');
  %set(get(h_c,'ylabel'),'String', 'n_b  [C/cm^3]');
  set(get(h_c,'ylabel'),'String', 'v_{EB} [m/s]');

  % add line plot on axis
  subplot(2,3,6);
%  subplot(2,1,2);
  plot(ZZ, qp(n_3D_counter).v_ExB(:,(2^INDX/2)+1));
  grid on;
  xlabel('z [um]')
  ylabel('E x B / |B|^2 [m/s]')
%  qp_FEX_max = max(qp_FEX_max, max(qp(n_3D_counter).v_ExB(:,(2^INDX/2)+1)));
%  qp_FEX_min = min(qp_FEX_min, min(qp(n_3D_counter).v_ExB(:,(2^INDX/2)+1)));
  myaxis = axis;
%  axis([min(ZZ)  max(ZZ)  qp_FEX_min  qp_FEX_max+eps]);
%  FEX_max_0 = max(max(qp(1).v_ExB(1:round(end*2/3), 1:end)));
  if(do_fixed_axes)
    axis([min(ZZ)  max(ZZ) -70 70]); % NJP
%    axis([min(ZZ)  max(ZZ) -50 50]);
%    axis([min(ZZ)  max(ZZ) -FEX_max_0 FEX_max_0]); % Dx-3
  end% if
  ExB_max = max(max(qp(n_3D_counter).v_ExB(1:round(end*2/3), 1:end)));
  ExB_min = min(min(qp(n_3D_counter).v_ExB(1:round(end*2/3), 1:end)));
%  title(['| ExB | = ' num2str(FEX_max, '%.1f') 'GV/m']);
  %title(['E_{x,max} = ' num2str(FEX_max, '%.1f') 'GV/m'  ', ' 'E_{x,min} = ' num2str(FEX_min, '%.1f') 'GV/m']);

%  pause;
end% if plot





if(do_FOCUSAXIS)
  subplot(1,1,1);
% find focusing field axis
x = size(qp(n_3D_counter).FEX, 2);
n_cut = Nx/2;
n_cut = 9;
%n_cut = 39;
FEX_zoom = abs(qp(n_3D_counter).FEX(:, (Nx/2-n_cut+1):(Nx/2+n_cut))');
[Y,I] = min(FEX_zoom, [],1);
I = I - n_cut - 1 + Nx/2; % adjust for new zero
%plot(I);
%pcolor(FEX_zoom);
%shading('flat');

% plot together with beam axis
plot(slice.z(n_beam, 1:end-1), slice.mean_x, '-ob');
hold on;
ZZ = (1:length(I)) * scale_z - offset_z0;
plot(ZZ, I * scale_x-offset_x0, '-xr');
hold off;
grid on;
xlabel('z [um]');
ylabel( plane_ylabel);
title(['s=' num2str(s_timestep*100, 3) 'cm, Step:' num2str(n_3D_timestep) 'DT', ', Slice:' qp_dump_plane ' ']);
legend('beam centroid', 'zero focus axis');
myaxis = axis;
axis([-110 0 -45 45 ]);
end% if  

% TEMP STOP ExB
%stop


%
%  plot beam phase space
%
if(do_beam_phase_space_plot )

% total beam analysis - merge beams (not stored in main qp var for space reasons)
qp_BEAMS = [];
for(n_beam=1:size(qp(n_3D_counter_beam).PP, 2))
  qp_BEAMS = [qp_BEAMS; qp(n_3D_counter_beam).PP(n_beam).BEAM];
end% for

% x, E histogram;  x : put 5 initial sigma
subplot(2,3,[1 4]);
% TEMP IPAC'12
%figure1 = figure('Color',[1 1 1]);
%subplot(1,1,1);
%set(0,'defaultaxesfontsize',24);


%hist_var1 = 3; % abscissa
hist_var1 = 1; % abscissa
hist_var2 = 6; % ordinate
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
E_min;
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
  E_min = 0; % FACET
  E_max = 40;% FACET
end% if
x_min = -100;
x_max= 100;
% IPAC'12
%x_min = -300;
%x_max= 300;
%  E_min = 0; % FACET
%  E_max = 40;% FACET
%  E_min = -10000; % x xp
%  E_max = 10000;% x xp
xedges = linspace(x_min, x_max, 512+0*round((x_max-x_min)*scale_x));
yedges = linspace(E_min, E_max, 512);
histmat = hist2(qp_BEAMS(:,hist_var1), qp_BEAMS(:,hist_var2), xedges, yedges);
h_g = pcolor(xedges,yedges,histmat); 
%h_c = colorbar; % ugly this big bar
if(do_fixed_axes)
  caxis([1 100]); % NJP
  caxis([0 10]);
  caxis([0 20]); % positrons
end% if
  caxis([0 5]);
shading('flat');
axis([x_min x_max E_min E_max+1e-10]); % for comp with AAC 2010 figure
%axis square tight;
xlabel(pplabel(hist_var1).var);
ylabel(pplabel(hist_var2).var);
%set(get(h_c,'ylabel'),'String', 'n_b  [a.u.]');
title(['s=' num2str(s_timestep*100, 3) 'cm, Step:' num2str(n_3D_timestep) 'DT', ', Slice:' qp_dump_plane ' ']);

%pause;

%  display(['Beam ' n_beam_str ': ' 'sigma_x=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_x(n_3D_counter_beam)) '  sigma_y=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_y(n_3D_counter_beam)) '  sigma_z=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_z(n_3D_counter_beam))  '  sigma_xp=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_xp(n_3D_counter_beam))  '  sigma_yp=' num2str(qp(n_3D_counter_beam).PP(n_beam).sigma_yp(n_3D_counter_beam))   '  mean_E=' num2str(qp(n_3D_counter_beam).PP(n_beam).mean_E(n_3D_counter_beam))     '  std_E / mean_E=' num2str(qp(n_3D_counter_beam).PP(n_beam).std_E(n_3D_counter_beam) / qp(n_3D_counter_beam).PP(n_beam).mean_E(n_3D_counter_beam)) ]);

end% if phase-space plot







%
%  plot Eprof
%
if(do_Eprof)
  n_hist = 101;
  ax_E_min = 1e10;
  ax_E_max = 0;
  n_max = 0;
  Q_max = 0;
  % extra loop to get n_max before plotting
  for nn_beam =1:size(qp(1).PP, 2),
    B = qp(n_3D_counter_beam).PP(nn_beam).BEAM;
    [n,E] = hist(B(:,6), n_hist);
    if( qp(1).PP(nn_beam).N0 > Q_max )
      Q_max = qp(1).PP(nn_beam).N0;
    end; % if
  end% for
  for nn_beam =1:size(qp(1).PP, 2),
    B = qp(n_3D_counter_beam).PP(nn_beam).BEAM;
    [n,E] = hist(B(:,6), n_hist);
    Q_frac = qp(1).PP(nn_beam).N0 / Q_max;
    if( max(n)*Q_frac > n_max )
      n_max = max(n)*Q_frac;
    end; % if
    if( max(n)*Q_frac > Eprof_n_max_glob )
      Eprof_n_max_glob = max(n)*Q_frac;
    end; % if
  end% for
  for nn_beam =1:size(qp(1).PP, 2),
  B = qp(n_3D_counter_beam).PP(nn_beam).BEAM;
  mean_E = mean(B(:,6));
  std_E = std(B(:,6));
  [n,E] = hist(B(:,6), n_hist);
  Q_frac = qp(1).PP(nn_beam).N0 / Q_max;
  E0 = qp(n_3D_counter).PP(nn_beam).gamma0*0.511e-03;
  if( min(E) < ax_E_min )
    ax_E_min = min(E);
  end; % if
  if( max(E) > ax_E_max )
    ax_E_max = max(E);
  end; % if
  n_acc_frac = n_acc_frac
  n_acc_T2_frac = n_acc_T2_frac
  n_dec_frac = n_dec_frac
    
  if(do_Eprof_allview)
    subplot(1,1,1);
  else
  subplot(2,3,3);
  end% if
   colorbar('off');  
%  hh = semilogy([mean_E mean_E+eps], [1 1e8], '-k');  
%  set(hh,'LineWidth',5)
 % hh = plot([mean_E+std_E mean_E+std_E+eps], [1 1e8], '-g');  
 % set(hh,'LineWidth',2)
 % hh = plot([mean_E-std_E mean_E-std_E+eps], [1 1e8], '-g');  
 % set(hh,'LineWidth',2)
%  hh = semilogy(E,n / max(n) * 1e4, ['-'
%  char(colorlist(nn_beam))]);
% pad with zeero-plot points
E = [min(E)-eps E max(E)+eps];
n = [0 n 0];
%  hh = plot(E,n /n_max * Q_frac * 1e3, ['-' char(colorlist(nn_beam))]);
  hh = plot(E,n /15722 * Q_frac * 1e3, ['-' char(colorlist(nn_beam))]);
  set(hh,'LineWidth',3)
  hold on;
end% for  
  hold off;
  grid on;
  xlabel('E [GeV]');
  ylabel('counts [a.u.] ');
  axis([ax_E_min ax_E_max*1.05 1 1e3*1.05]);
  % 
  set(0,'defaulttextfontsize',18);
%  set(0,'defaulttextfontname','Courier');
%  text(31, 100, ['Li-%: ' sprintf('%0.2g', n_dec_frac*100)]);
%  hh = legend(['<E>=' num2str(mean_E, 2) 'GeV'], ['\sigma_E/E=' num2str(std_E/mean_E*100, 2) '%']);
%  set(hh,'FontSize',12);
  if(do_Eprof_allview)
    ht = title(['   s=' sprintf('%.0f', qp(n_3D_counter).s_timestep*1e2) ' cm' ',   <E>_{WB}=' sprintf('%.1f', qp(end).PP(2).mean_E) ' GeV' ', \sigma /E_{WB}=' sprintf('%.1f', qp(end).PP(2).std_E/qp(end).PP(2).mean_E*100) ' %' ',   <E>_{DB}=' sprintf('%.1f', qp(end).PP(1).mean_E) ' GeV   ']);
    axis([0 11 0 1050]);
  pos=get(ht,'Position');
  pos(1)=pos(1)+0.3;
  pos(2)=pos(2)-40;
set(ht,'Position',pos);
  else
    title(['s = ' num2str(qp(n_3D_counter).s_timestep*1e2, 3) ' cm']);
  end% if
%  set(hh,'FontName','Courier');
if(size(qp(1).PP, 2) == 2)
  legend('DB', 'WB', 'Location', 'North');
  myaxis = axis;
  axis([0 myaxis(2) myaxis(3) myaxis(4)]);
end% if
end% if



if(~do_movie && do_both_plot_types)
  pause;
end% if

% movie-maker
if( do_movie )
for(n_dt_dt_per_frame=1:N_dt_per_frame)
  n_movframe = n_movframe + 1;
  mymovdir = [working_dir 'movie_frames/'];
  myfile_mov = [mymovdir 'frame' num2str(n_movframe, '%5.5d') '.png'];
  saveas(gcf, myfile_mov, 'png');
end% for
end% if



end% 3D timestep loop

end% if plot


