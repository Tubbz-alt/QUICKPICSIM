function my_gen_rpinput(myfilein, myfileout, Plasma_Density, plasma_Z, plasma_PREION, plasma_s_prop, beam_charge, beam_mean_E, beam_sigma_E_E, beam_sigma_x, beam_sigma_y, beam_sigma_z, beam_emitt_x, beam_emitt_y, BEAM_EV, beam_match, emitt_match, beam_z_pos, tilt_x, tilt_y, histZ, histZcount)

working_dir = '/Users/eadli/Dropbox/SLAC/quickpic/QUICKPICSIM/'; eval(['run ' working_dir 'my_SI_params.m']); % import my standard SI constants
fun_dir = '/Users/eadli/Dropbox/SLAC/quickpic'; addpath(fun_dir);

if( nargin < 17 )
  error('EA ERR: Not enough input arguments to "my_gen_rpinput".');
end% if
if( nargin < 18 )
  beam_z_pos = -1; % put bunch(es) at Box_Z/2
end% if
if( nargin < 19 )
  tilt_x = [0 0 0]';
end% if
if( nargin < 20 )
  tilt_y = [0 0 0]';
end% if
beam_z_is_custom = 1;
if( nargin < 22 )
  beam_z_is_custom = 0;
end% if


%
% number of beams: length of beam parameter vectors
%
N_beams = size(beam_charge, 2)


%
% derived and auxiliary input for quickpic
%
N_PHASESPACE = 17; % 2^N particles in phase space output
% these are "dangerous" too touch ... simulations may crash
NPX =  7; % 2^N
NPY = 7; % 2^N
NPZ = 8; % 2^N
%
omega_p = sqrt(Plasma_Density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi % plasma wavelength
k_p = 2*pi/lambda_p;
%
% Box size: max( 2.5 times Bubble Radius, 10 sigma_x, 10 sigma_y)
gamma = beam_mean_E/1e6 / 0.511
if( beam_match )
  beam_sigma_x = sqrt(beam_emitt_x ./ k_p .*sqrt(2./gamma) )
  beam_sigma_y = sqrt(beam_emitt_y ./ k_p .*sqrt(2./gamma) )
end% if
if( emitt_match )
  k_p
  beam_sigma_x
  gamma
  beam_emitt_x = k_p*beam_sigma_x.^2*sqrt(gamma./2)
  beam_emitt_y = k_p*beam_sigma_y.^2*sqrt(gamma./2)
end% if
if(beam_z_is_custom)
  % calculate sigma_z from distribution 
  for(nn_bunch = 1:size(histZ, 1))
    beam_sigma_z(nn_bunch) = my_std_lp(histZ(nn_bunch, :), histZcount(nn_bunch, :))
    beam_sigma_z = 500e-6;
    % well, but for peak current, use rather peak current directly!
    dz = histZ(nn_bunch, 2) - histZ(nn_bunch, 1)
    lambda_peak = max(histZcount(nn_bunch, :)) / sum(histZcount(nn_bunch, :))/dz
  end% for
  % well, but for peak current, use rather peak dist directly!
  I_peak = beam_charge*SI_c * lambda_peak
else
  I_peak = beam_charge*SI_c ./ sqrt(2*pi*beam_sigma_z.^2) % Gaussian
end% if
R_bubble = (1/0.84) * 2 / k_p * sqrt(max(I_peak) / 17e3)
%Box_X_exact = max(1.0 * 4 * R_bubble, max(10*max(beam_sigma_x), 10*max(beam_sigma_y))) % [m]
Box_X_exact = 1.8*max(1.0 * 4 * R_bubble, max(5*max(beam_sigma_x), 5*max(beam_sigma_y))) % [m]
% if tilted beams add space in x
if( (sum(abs(tilt_x))+sum(abs(tilt_y))) > 0)
  Box_X_exact = Box_X_exact +2*sqrt(3)*beam_sigma_z*max( max(abs(tilt_x(2,:))), max(abs(tilt_y(2,:))));
end% if
Box_X = 1.5 * ceil(Box_X_exact*1e6/10)*10 % [um] - round the box size
Box_X = max(Box_X); % if multiple bunches, use largest Box_X calc
% make odd number of transverse cells (in order to center the beam)
Box_X = Box_X + 1; % if multiple bunches, use largest Box_X calc
%Box_X = 580
Box_Y = Box_X % [um]
Box_Z_exact = max(6.0*max(beam_sigma_z), 2.5 * lambda_p) % [um]
Box_Z = round(Box_Z_exact*1e6/20)*20 % [um] - round the box size
%Box_Z = 190
Box_Z = Box_Z + 1;
%
% TEMP BOX REDUCE
%
%Box_X = Box_X * 0.6;
%Box_Y = Box_Y * 0.6;
%Box_Z = Box_Z * 0.6;
skin_frac = 1/20; % grid size in fractions of skin-depth
INDX = floor( log(Box_X*1e-6 * k_p / skin_frac) / log(2) )
INDY = floor( log(Box_Y*1e-6 * k_p / skin_frac) / log(2) )
INDZ = floor( log(Box_Z*1e-6 * k_p / skin_frac) / log(2) )
%INDX = 9 % temp
%INDY = 9 % temp
%Box_Z = Box_Z - 50 % temp
% force min 2^6 per dim
if( INDX < 6 )
  INDX = 6;
end% if
if( INDY < 6 )
  INDY = 6;
end% if
if( INDZ < 6 )
  INDZ = 6;
end% if
% enfore cellZ < cellX (other effects have not been predicted) 
while( (Box_Z/2^INDZ) > (Box_X/2^INDX) )
  INDZ = INDZ + 1
end% while

% Where to put bunches?
%   ( if custom z then N/A )
if(~beam_z_is_custom)
  if( beam_z_pos == -1)
    % if user gives nothing, or -1, put (all..) bunches on Box_Z * 1/3
    beam_z_pos = round(Box_Z * 1/3*1.4) * ones(N_beams, 1)
  end% if
end% if
if( beam_z_pos == -1)
  beam_z_pos = round(Box_Z* 1/3*1.4) * ones(N_beams, 1);
end% if

% Where to put bunches transversally if tilt?
%   force bunch towards center of box
comp_tilt_x = - round(Box_Z/2 * tilt_x(2, :));
comp_tilt_y = - round(Box_Z/2 * tilt_y(2, :));


%
DT =round(sqrt(2*min(gamma))/10 / 1.5) % 1.5 is deceleration factor
if( BEAM_EV )
  TEND = floor(plasma_s_prop / (SI_c / omega_p))+0.1;
  DT_OUTPUT = 10; % output data every n'th timestep
  DT_RESTART = DT_OUTPUT*4;
else
  TEND = DT+0.1;
  DT_OUTPUT = 1; % output data every timestep
  DT_RESTART = 0;
end% if
Plasma_Density_str = num2str(round(Plasma_Density), '%.3G');
Plasma_Density_str = strrep(Plasma_Density_str,'E+0','E');
Plasma_Density_str = strrep(Plasma_Density_str,'E+','E');
TEND_str = num2str((TEND), '%10G');
TEND_str = strrep(TEND_str,'E+0','E');
TEND_str = strrep(TEND_str,'E+','E')
%
store_QEB_3D = 0; % can give total charge, but very resource and time consuming

% stages: here assume n=128 jobs (must be specified at job set-up)
if( BEAM_EV )
  n_tasks = 128;
else
  n_tasks = 8; % use this for express jobs
end% if
Num_Stages = max(n_tasks / 2^(INDX-3), 1);



%  conversions
% mean energy  [cannot import energy spread]
% hist Z
% x, emitx, y, emit
% tilt in x
% tilt in y
qp_emitt_x = beam_emitt_x*1e6
qp_emitt_y = beam_emitt_y*1e6
qp_sigma_x = beam_sigma_x*1e6
qp_sigma_y = beam_sigma_y*1e6
qp_e_spread = beam_sigma_E_E; 
qp_e_spread = zeros(1,N_beams); % Feature not working, so set to 0
qp_num_part = beam_charge / SI_e;

% custom z distribution
if(~beam_z_is_custom)
  qp_sigma_z = beam_sigma_z*1e6;
end% if


%
% 
%

%
% generate rpinput file by replacing modified lines in a standard rpinput
%
fid = fopen(myfilein, 'r');
fidout = fopen(myfileout, 'w');
section_beam = 0;
section_plasma = 0;
section_neutral = 0;
section_species = 0;
section_simsys = 0;
section_nbeams = 0;
section_simtime = 0;
section_fielddiag = 0;
section_beamdiag = 0;
section_plasmadiag = 0;
section_beamphasediag = 0;
section_restart = 0;
row = 0;
n_row = 0;
while( row ~= -1 )
  n_row = n_row + 1;
  row = fgets(fid);
  if( row ~= -1 )
    % TAG SECTIONS
    if( ~isempty(strfind(row, '&Pipeline')) )
      section_pipeline = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_pipeline = 0;
    end% if
    if( ~isempty(strfind(row, '&Simulation_Sys')) )
      section_simsys = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_simsys = 0;
    end% if
    if( ~isempty(strfind(row, '&Num_Beams')) )
      section_nbeams = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_nbeams = 0;
    end% if
    if( ~isempty(strfind(row, '&Beam')) && isempty(strfind(row, '&Beam_')) )
      section_beam = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_beam = 0;
    end% if
    if( ~isempty(strfind(row, '&Plasma')) && isempty(strfind(row, '&Plasma_')) )
      section_plasma = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_plasma = 0;
    end% if
    if( ~isempty(strfind(row, '&Neutral')) && isempty(strfind(row, '&Neutral_')) )
      section_neutral = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_neutral = 0;
    end% if
    if( ~isempty(strfind(row, '&Species')) )
      section_species = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_species = 0;
    end% if
    if( ~isempty(strfind(row, '&Simulation_time')) )
      section_simtime = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_simtime = 0;
    end% if
    if( ~isempty(strfind(row, '&Field_Diag')) )
      section_fielddiag = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_fielddiag = 0;
    end% if
    if( ~isempty(strfind(row, '&Beam_Diag')) )
      section_beamdiag = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_beamdiag = 0;
    end% if
    if( ~isempty(strfind(row, '&Plasma_Diag')) )
      section_plasmadiag = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_plasmadiag = 0;
    end% if
    if( ~isempty(strfind(row, '&Beam_Phase_Space_Diag')) )
      section_beamphasediag = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_beamphasediag = 0;
    end% if
    if( ~isempty(strfind(row, '&Restart_File')) )
      section_restart = 1;
    end% if
    if( ~isempty(strfind(row, '/')) )
      section_restart = 0;
    end% if
    

       % STAGES SECTION
   if( ~isempty(strfind(row, 'Num_Stages') ) && section_pipeline == 1)
      fprintf(fidout, [' Num_Stages = ' num2str(Num_Stages) '\n']  );
      
       % SIMSYS SECTION
   elseif( ~isempty(strfind(row, 'Box_X') ) && section_simsys == 1)
      fprintf(fidout, [' Box_X=' num2str(Box_X) ', Box_Y=' num2str(Box_Y) ', Box_Z=' num2str(Box_Z) ',' '\n']  );
    elseif( ~isempty(strfind(row, 'INDX') ) && section_simsys == 1)
      fprintf(fidout, [' INDX =  ' num2str(INDX) ' , INDY = ' num2str(INDY) ', INDZ = ' num2str(INDZ) '\n']  );

       % NUMBEAMS SECTION
   elseif( ~isempty(strfind(row, 'NBeams') ) && section_nbeams == 1)
      fprintf(fidout, [' NBeams = ' num2str(N_beams) '\n']  );
      
      % BEAM SECTION(S)
   elseif( section_beam == 1)
for(n_beam=1:N_beams),
fprintf(fidout, ['\n']  );
fprintf(fidout, ['&Beam' '\n']  );
fprintf(fidout, [' BEAM_EVOLUTION = .true.' '\n']  );
if( beam_z_is_custom ),
  % for user z-distribution, this number must be "large"
  fprintf(fidout, [' MIN_BEAM_PARTICLE = 5000000' '\n']  );
else
%  fprintf(fidout, [' MIN_BEAM_PARTICLE = 8' '\n']  );
   fprintf(fidout, [' MIN_BEAM_PARTICLE = 5000000' '\n']  ); % update after input from W. An 2012-05-14
end% if
fprintf(fidout, [' NPX =  ' num2str(2^NPX) ', NPY = ' num2str(2^NPY) ', NPZ = ' num2str(2^NPZ) '\n']);
fprintf(fidout, [' Charge = -1.0' '\n']  );
fprintf(fidout, [' Mass = 1.0' '\n']  );
fprintf(fidout, [' Gamma = ' num2str(round(gamma(n_beam)), '%d'), ',' '\n']  );
qp_num_part_str = num2str(round(qp_num_part(n_beam)), '%.3G');
qp_num_part_str = strrep(qp_num_part_str,'E+0','E');
qp_num_part_str = strrep(qp_num_part_str,'E+','E')
fprintf(fidout, [' Num_Particle = ' qp_num_part_str ',' '\n']);
fprintf(fidout, [' VDX =   0.0, VDY =   0.0, VDZ =  0.0' '\n']  );
if( beam_z_is_custom )
  fprintf(fidout, [' Init_Routine = 3' '\n']  );
else
  fprintf(fidout, [' Init_Routine = 1' '\n']  );
end% if
fprintf(fidout, [' BEAM_PROFILE = ''test.hdf''' '\n']  );
fprintf(fidout, [' QUIET_START = .true.' '\n']  );
fprintf(fidout, [' Parameter_Array(1:1,1:3) = ' sprintf('%.3f', Box_X/2 + 0*comp_tilt_x(n_beam)) ',' sprintf('%.3f', Box_Y/2 + comp_tilt_y(n_beam))  ',' num2str(beam_z_pos(n_beam), '%.2f') '\n']);
if( beam_z_is_custom )
  fprintf(fidout, [' Parameter_Array(2:2,1:3) = ' num2str(qp_sigma_x(n_beam), '%.2f') ',' num2str(qp_sigma_y(n_beam), '%.2f')  ',' num2str(length(histZcount), '%.0d') '\n']);
else
  fprintf(fidout, [' Parameter_Array(2:2,1:3) = ' num2str(qp_sigma_x(n_beam), '%.2f') ',' num2str(qp_sigma_y(n_beam), '%.2f')  ',' num2str(qp_sigma_z(n_beam), '%.2f') '\n']);
end% if
if( beam_z_is_custom)
  fprintf(fidout, [' Parameter_Array(3:3,1:9) = ' num2str(qp_emitt_x(n_beam), '%.3f') ',' num2str(qp_emitt_y(n_beam), '%.3f') ',' num2str(qp_e_spread(n_beam), '%.3f') ',' num2str(tilt_x(3, n_beam), '%.3f') ',' num2str(-tilt_x(2, n_beam), '%.3f') ',' num2str(tilt_x(1, n_beam), '%.3f') ',' num2str(tilt_y(3, n_beam), '%.3f') ',' num2str(-tilt_y(2, n_beam), '%.3f') ',' num2str(tilt_y(1, n_beam), '%.3f')  '\n'] );
else
  fprintf(fidout, [' Parameter_Array(3:3,1:3) = ' num2str(qp_emitt_x(n_beam), '%.3f') ',' num2str(qp_emitt_y(n_beam), '%.3f') ',' num2str(qp_e_spread(n_beam), '%.3f')  '\n'] );
end% if
if( beam_z_is_custom)
  histZcount_str = [];
  for ( n=1:length(histZcount)-1),
    histZcount_str = [histZcount_str num2str(histZcount(n_beam, n), '%.1f') ','];
  end% for
  histZcount_str = [histZcount_str num2str(histZcount(n_beam, n+1), '%.1f')];
  histZ_str = [];
  for ( n=1:length(histZ)-1),
    histZ_str = [histZ_str num2str(histZ(n_beam, n)*1e6 + round(Box_Z/2)-0, '%.1f') ',']; %  - sqrt(3)*beam_sigma_z(n_beam)*1e6
  end% for
  histZ_str = [histZ_str num2str(histZ(n_beam, n+1)*1e6 + round(Box_Z/2)-0, '%.1f')];
  fprintf(fidout, [' Parameter_Array(4:4,1:' num2str(length(histZcount)) ') = ' histZcount_str  '\n'] );
else
  fprintf(fidout, [' Parameter_Array(4:4,1:3) = '  num2str(tilt_x(3, n_beam), '%.3f') ',' num2str(-tilt_x(2, n_beam), '%.3f') ',' num2str(tilt_x(1, n_beam), '%.3f')  '\n'] );
end% if
if( beam_z_is_custom)
  fprintf(fidout, [' Parameter_Array(5:5,1:' num2str(length(histZ)) ') = ' histZ_str  '\n'] );
else
  fprintf(fidout, [' Parameter_Array(5:5,1:3) = '  num2str(tilt_y(3, n_beam), '%.3f') ',' num2str(-tilt_y(2, n_beam), '%.3f') ',' num2str(tilt_y(1, n_beam), '%.3f')  '\n'] );
end% if
fprintf(fidout, [' Use_Shifter = .false.' '\n']  );
fprintf(fidout, [' Shifter_Nsec = 4' '\n']  );
fprintf(fidout, [' Shifter_Parameter(1:1,1:4) = 0.,0.,1.5,0.' '\n']  );
fprintf(fidout, [' Shifter_Parameter(2:2,1:4) = 0.,0.,0.,0.' '\n']  );
fprintf(fidout, [' Shifter_Parameter(3:3,1:4) = 0.,78.,155.1,155.2' '\n']  );
fprintf(fidout, [' Use_Destroyer = .true.' '\n']  );
fprintf(fidout, [' Destroyer_NCriteria = 5' '\n']  );
fprintf(fidout, [' Destroyer_Criteria(1:1,1:5)=1,1,2,2,6' '\n']  );
fprintf(fidout, [' Destroyer_Criteria(2:2,1:5)=0,' num2str(Box_X-2, '%.0f') ',0,' num2str(Box_Y-2, '%.0f')  ',0' '\n']);
fprintf(fidout, [' Destroyer_Criteria(3:3,1:5)=2,' num2str(Box_X, '%.0f') ',2,' num2str(Box_Y, '%.0f')  ',100' '\n']);
fprintf(fidout, [' Use_Radiation_Damping = .false.' '\n']  );
fprintf(fidout, ['/' '\n']  );
end% for
section_beam = 0;


     % PLASMA SECTION
   elseif( ~isempty(strfind(row, 'Plasma_Density') ) && section_plasma == 1)
      fprintf(fidout, [' Plasma_Density=' Plasma_Density_str '\n']  );
    elseif( ~isempty(strfind(row, 'Plasma_Density') ) && section_plasma == 1)
      fprintf(fidout, [' Plasma_Density=' Plasma_Density_str '\n']  );
   elseif( ~isempty(strfind(row, 'Nspecies') ) && section_plasma == 1)
     if( plasma_PREION == 1)
       fprintf(fidout, [' Nspecies=1' '\n']  );
     else
       fprintf(fidout, [' Nspecies=0' '\n']  );
     end% if
   elseif( ~isempty(strfind(row, 'Nneutral') ) && section_plasma == 1)
     if( plasma_PREION == 1)
       fprintf(fidout, [' Nneutrals=0' '\n']  );
     else
       fprintf(fidout, [' Nneutrals=1' '\n']  );
     end% if

     % NEUTRAL SECTION
    elseif( ~isempty(strfind(row, 'Neutral_gas') ) && section_neutral == 1)
      fprintf(fidout, [' Neutral_gas = ' num2str(plasma_Z) '\n']  );

     % SPECIES SECTION 
   elseif( ~isempty(strfind(row, 'NP2') ) && section_species == 1)
     n_plasma_particles_per_cell = 4; % form UCLA
     fprintf(fidout, [' NP2 = ' num2str(round(sqrt((2^INDX)^2 * n_plasma_particles_per_cell))) '\n']  );

      % SIMTIME SECTION
    elseif( ~isempty(strfind(row, 'TEND') ) && section_simtime == 1) 
      fprintf(fidout, [' TEND =' TEND_str ', DT = ' num2str(DT, '%.1f') '  ,' '\n']  );

      % DIAG SECTION
    elseif( ~isempty(strfind(row, 'DFESLICE') ) && section_fielddiag == 1) 
      fprintf(fidout, [' DFESLICE=' num2str(DT_OUTPUT) ', EX0=0, EY0=' num2str(round(Box_Y/2), '%.0f') ', EZ0=0' '\n']  );
    elseif( ~isempty(strfind(row, 'DFBSLICE') ) && section_fielddiag == 1) 
      fprintf(fidout, [' DFBSLICE=' num2str(DT_OUTPUT) ', BX0=0, BY0=' num2str(round(Box_Y/2), '%.0f') ', BZ0=0' '\n']  );
   elseif( ~isempty(strfind(row, 'DFQEB')) && isempty(strfind(row, 'DFQEBSLICE'))  && section_beamdiag == 1)
     if(store_QEB_3D)
       fprintf(fidout, [' DFQEB=' num2str(DT_OUTPUT) ',' '\n']  );
     else
       fprintf(fidout, [' DFQEB=0,' '\n']  );
     end% if
    elseif( ~isempty(strfind(row, 'DFQEBSLICE'))  && section_beamdiag == 1)
      fprintf(fidout, [' DFQEBSLICE=' num2str(DT_OUTPUT) ' , QEBX0=0., QEBY0=' num2str(round(Box_Y/2), '%.0f') ',  QEBZ0=0' '\n']  );
    elseif( ~isempty(strfind(row, 'DFQEPSLICE'))  && section_plasmadiag == 1)
      fprintf(fidout, [' DFQEPSLICE=' num2str(DT_OUTPUT) ' , QEPX0= 0, QEPY0=' num2str(round(Box_Y/2), '%.0f') ',  QEPZ0=0' '\n']  );
    elseif( ~isempty(strfind(row, 'DUMP_PHA_BEAM'))  && section_beamphasediag == 1)
      fprintf(fidout, [' DUMP_PHA_BEAM=.true., DFPHA_BEAM=' num2str(DT_OUTPUT) ',' '\n']  );
    elseif( ~isempty(strfind(row, 'DSAMPLE_BEAM'))  && section_beamphasediag == 1)
      fprintf(fidout, [' DSAMPLE_BEAM = ' num2str(2^(NPX+NPY+NPZ-N_PHASESPACE)) '\n']  );

      % RESTART SECTION
   elseif( ~isempty(strfind(row, 'DUMP_RST_FILE'))  && section_restart == 1)
     if( DT_RESTART)
       fprintf(fidout, [' DUMP_RST_FILE = .true.,  DFRST=' num2str(DT_RESTART) '\n']  );
     else
       fprintf(fidout, [' DUMP_RST_FILE = .false.,  DFRST=1' '\n']  );
     end% if
   else
      fprintf(fidout, row);
    end%if
  end%if
end% while
fclose(fid);
fclose(fidout);
