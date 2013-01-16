function [b_mat, a_mat, b_0, a_0, b_min, a_min, b_max, a_max, b_beat_rel, z, beta_evol] =my_calc_ramp(n0, gamma, d_beta_beta0, do_plot, z_ramp)


%
% Plasma ramp stuff
%
%qmat = csvread('/Users/eadli/Dropbox/SLAC/elegant/PLASMALENS/ElegantPlasma/qmat.csv');  
% now: instead generate your plasma ramp from scratch

working_dir = '/Users/eadli/Dropbox/SLAC/quickpic/QUICKPICSIM/'; eval(['run ' working_dir 'my_SI_params.m']); % import my standard SI constants> 

omega_p = sqrt(n0*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
k_p_inv = 1/ k_p; % skin-depth and UCLA length scale

% input beta
beta_1 = sqrt(2*gamma) / k_p;

%%% plasma ramp parameters
roll_up = 1; % 0: roll_down

% speed of plasma rollup
w      = 5 / z_ramp; % [m^-1]
%w      = 33;
nsteps = 2048;
%nsteps = 64;
%nsteps = 128;
%nsteps = 47; % for quickpic two ramp
% step size
dz     = z_ramp/nsteps; % [m]

z = dz*(0:(nsteps-1));
%np = n0*1e6*(tanh(w*(z-z_ramp/2))+1);
np =  n0*1e6*(tanh(w*(z-z_ramp/2))+1)/2; % corrected (was factor 2 too large)
%np = n0*(1-tanh(w (z-z0/2)) ) ./ (1- tanh(w*(-z0/2)) ); % JF's version of fit
%plot(np)
%stop
%min(np)%max(np)
%max(np)



%
% manually scale to match 0 to 1 (fit function is cut before 0 and 1)
np = np - np(1)+eps;
scale_f = (n0*1e6)/np(end);
np = np * scale_f;


%
% write ramp for input into quickpic
%   [s [um]    np/n0 [-]
do_write_np = 1;
flat_top_s = 0.01; % [m]
flat_top_s = 0.30; % [m]
if(do_write_np)
  filename = '/tmp/lens_np.txt';
  fid = fopen(filename, 'w');
  Nsec = nsteps+3;
  if(flat_top_s > 0.0)
    Nsec = Nsec + 1;
  end% if
 fprintf(fid, ' Density_Variation=.true.\n');
 fprintf(fid, ' Density_Variation_NSec=%d\n', Nsec);
 fprintf(fid, ' Density_Variation_Fs(1:%d) = ', Nsec);
  % define initial flat top
    if(roll_up)
  if(flat_top_s > 0.0)
    fprintf(fid, '%.3e,', 1.0e-10);
  end% if
    else
  fprintf(fid, '%.3e,', 1.0);
  if(flat_top_s > 0.0)
    fprintf(fid, '%.3e,', 1.0);
  end% if
      end% if
  for n=1:nsteps
    if(roll_up)
      n_pick = n;
    else
      n_pick = nsteps-n+1;
    end;
    myvalue = np(n_pick) / n0 / 1e6;
    if(myvalue < 1e-10)
      myvalue = 1e-10; % "light paranoia", avoiding zero values for QuickPIC
    end% if
    fprintf(fid, '%.3e,', myvalue);
  end% if
  % define final points of low density (5 cm of ~ zero dens see drift)
    if(roll_up)
  fprintf(fid, '%.3e,', 1.0);
  fprintf(fid, '%.3e', 1.0);
    else
  fprintf(fid, '%.3e,', 1.0E-10);
  fprintf(fid, '%.3e', 1.0E-10);
      end% if
    fprintf(fid, '\n');
  % define initial flat top
  fprintf(fid, ' Density_Variation_s(1:%d) = ', Nsec);
  fprintf(fid, '%.3e,', 0.0);
  if(flat_top_s > 0.0)
    if(roll_up)
      % do nothing
    else
          fprintf(fid, '%.3e,', flat_top_s*1e6);
    end% if
  end% if
  for n=1:nsteps
    if(roll_up)
      n_pick = n;
    else
      n_pick = nsteps-n+1;
    end;
    L = z_ramp/nsteps;
    if(roll_up)
       roll_up_delta = 50e3;
      fprintf(fid, '%.3e,', L*n*1e6);
    else
      fprintf(fid, '%.3e,', L*n*1e6 + flat_top_s*1e6);
    end% if
  end% if
  % define final points of low density (5 cm of ~ zero dens see drift)
    if(roll_up)
      fprintf(fid, '%.3e,', L*n*1e6 + L*1e6);
      fprintf(fid, '%.3e',  L*n*1e6 + L*1e6 + flat_top_s*1e6);
    else
      fprintf(fid, '%.3e,', L*n*1e6 + flat_top_s*1e6 + L*1e6);
      fprintf(fid, '%.3e', L*n*1e6 + flat_top_s*1e6 + L*1e6 + 50e3);
    end% if
    fprintf(fid, '\n');
  fclose(fid);
end% if  


% classical radius of elec.
re     = SI_re; % [m]
K_focus = ( sqrt(2*pi*SI_re*np/gamma) );
% K_focus = k_p / sqrt(2*gamma); % equiv to above
qmat = [];
for n=1:nsteps
  qmat(n,:) = [1 1 cos(K_focus(n)*dz) 1 sin(K_focus(n)*dz)/K_focus(n)  1   -K_focus(n)*sin(K_focus(n)*dz) 1 cos(K_focus(n)*dz)];
end

do_write_elegant = 1;
if(do_write_elegant)
  L = z_ramp/nsteps;
  filename = '/tmp/lens_elegant.par.txt';
  fid = fopen(filename, 'w');
  for n=1:nsteps
    if(roll_up)
      n_pick = n;
    else
      n_pick = nsteps-n+1;
    end;
    r11 = qmat(n_pick, 3);
    r12 = qmat(n_pick, 5);
    r21 = qmat(n_pick, 7);
    r22 = qmat(n_pick, 9);
    fprintf(fid, 'PLASMALENS %d L %.10e\n', n, L);
    fprintf(fid, 'PLASMALENS %d R11 %.10e\n', n, r11);
    fprintf(fid, 'PLASMALENS %d R33 %.10e\n', n, r11);
    fprintf(fid, 'PLASMALENS %d R12 %.10e\n', n, r12);
    fprintf(fid, 'PLASMALENS %d R34 %.10e\n', n, r12);
    fprintf(fid, 'PLASMALENS %d R21 %.10e\n', n, r21);
    fprintf(fid, 'PLASMALENS %d R43 %.10e\n', n, r21);
    fprintf(fid, 'PLASMALENS %d R22 %.10e\n', n, r22);
    fprintf(fid, 'PLASMALENS %d R44 %.10e\n', n, r22);
  end
  fclose(fid);
end% if  



do_write_elegant = 0;
if(do_write_elegant)
  L = z_ramp/nsteps;
  filename = '/tmp/lens_elegant.par.txt';
  fid = fopen(filename, 'w');
  for n=1:nsteps
    if(roll_up)
      n_pick = n;
    else
      n_pick = nsteps-n+1;
    end;
    r11 = qmat(n_pick, 3);
    r12 = qmat(n_pick, 5);
    r21 = qmat(n_pick, 7);
    r22 = qmat(n_pick, 9);
    fprintf(fid, 'LENSDRIFT1 %d L %.10e\n', n, L/2);
    f = -1 / r21;
    fprintf(fid, 'LENSTHIN %d FX %.10e\n', n, f);
    fprintf(fid, 'LENSTHIN %d FY %.10e\n', n, f);
    fprintf(fid, 'LENSDRIFT2 %d L %.10e\n', n, L/2);
  end
  fclose(fid);

  % drift instead of lens
  L = z_ramp/nsteps;
  filename = '/tmp/drift_elegant.par.txt';
  fid = fopen(filename, 'w');
  for n=1:nsteps
    if(roll_up)
      n_pick = n;
    else
      n_pick = nsteps-n+1;
    end;
    r11 = qmat(n_pick, 3);
    r12 = qmat(n_pick, 5);
    r21 = qmat(n_pick, 7);
    r22 = qmat(n_pick, 9);
    fprintf(fid, 'LENSDRIFT1 %d L %.10e\n', n, L/2);
    f = -1 / r21;
    fprintf(fid, 'LENSTHIN %d FX %.10e\n', n, f*1e10);
    fprintf(fid, 'LENSTHIN %d FY %.10e\n', n, f*1e10);
    fprintf(fid, 'LENSDRIFT2 %d L %.10e\n', n, L/2);
  end
  fclose(fid);
end% if  



do_write_elegant = 1;
if(do_write_elegant)
  L = z_ramp/nsteps;
  filename = '/tmp/lens_elegant_X.par.txt';
  fid = fopen(filename, 'w');
  for n=1:nsteps
    if(roll_up)
      n_pick = n;
    else
      n_pick = nsteps-n+1;
    end;
    K = K_focus(n_pick)^2;
    fprintf(fid, 'QUADLENS %d L %.10e\n', n, L);
    fprintf(fid, 'QUADLENS %d K1 %.10e\n', n, K);
  end
  fclose(fid);
  
  filename = '/tmp/lens_elegant_Y.par.txt';
  fid = fopen(filename, 'w');
  for n=1:nsteps
    if(roll_up)
      n_pick = n;
    else
      n_pick = nsteps-n+1;
    end;
    K = K_focus(n_pick)^2;
    fprintf(fid, 'QUADLENS %d L %.10e\n', n, L);
    fprintf(fid, 'QUADLENS %d K1 %.10e\n', n, -K);
  end
  fclose(fid);
end% if


K_focus_flattop = ( sqrt(2*pi*SI_re*n0*1e6/gamma) );
%K_focus_flattop = k_p / sqrt(2*gamma)  % equiv to above

qmat_flattop = [1 1 cos(K_focus_flattop*dz) 1 sin(K_focus_flattop*dz)/K_focus_flattop  1   -K_focus_flattop*sin(K_focus_flattop*dz) 1 cos(K_focus_flattop*dz)];

R_tot = eye(2,2);
for n=1:size(qmat,1),
  r11 = qmat(n, 3);
  r12 = qmat(n, 5);
  r21 = qmat(n, 7);
  r22 = qmat(n, 9);
  %
  R = [r11 r12; r21 r22];
  R_tot = R*R_tot;
end% for

alpha_1 = 0;
gamma_1 = (1+alpha_1)^2 / beta_1;

B1 = [beta_1 -alpha_1; -alpha_1 gamma_1];
B0 = inv(R_tot)*B1*inv(R_tot');

% add eventual error in beta to study beta beat
beta_0 = B0(1,1) * (1 + d_beta_beta0) * 1;
alpha_0 = -B0(1,2);
gamma_0 = (1+alpha_0^2) / beta_0;

% regeberate B0 with error in incoming beta
B0 = [beta_0 -alpha_0; -alpha_0 gamma_0];

B1_test = R_tot*B0*R_tot';

% backtrack; 
% beta evolution all the way (to find minimum)
%
B0_evol(:,:,1) = B0;
R_tot = eye(2,2);
for n=1:size(qmat,1),
  r11 = qmat(n, 3);
  r12 = qmat(n, 5);
  r21 = qmat(n, 7);
  r22 = qmat(n, 9);
  %
  R = [r11 r12; r21 r22];
  B0_evol(:,:,n+1) = R*B0_evol(:,:,n)*R';
end% for

b_mat = (B0_evol(1,1,end));
a_mat = ((B0_evol(1,2,end)));
b_0 = (B0_evol(1,1,1));
a_0 = ((B0_evol(1,2,1)));
b_min = min(B0_evol(1,1,:));
a_min = min((B0_evol(1,2,:)));
b_max = max(B0_evol(1,1,:));
a_max = max((B0_evol(1,2,:)));
if( roll_up )
  a_0 = -a_0;
end% if


% 
% continue beta evolution past ramp (to find beta beating)
%
%R_tot = eye(2,2);
% B0_evol(1,1,n+1)
% -B0_evol(2,1,n+1)
%B_mat_new = sqrt(2*gamma) / k_p 
%B_mat_new = (1 / K_focus_flattop)
%B0_evol(:,:,size(qmat,1)+1) = [B_mat_new 0; 0 1/B_mat_new];

% is difference of using real flat top significant?
%%% display
  r11 = qmat(end, 3);
  r12 = qmat(end, 5);
  r21 = qmat(end, 7);
  r22 = qmat(end, 9);
  r11 = qmat_flattop(3);
  r12 = qmat_flattop(5);
  r21 = qmat_flattop(7);
  r22 = qmat_flattop(9);
  R = [r11 r12; r21 r22];
  %
for n=(size(qmat,1)+1):(2*size(qmat,1)),
  r11 = qmat(end, 3);
  r12 = qmat(end, 5);
  r21 = qmat(end, 7);
  r22 = qmat(end, 9);
%  r11 = qmat_flattop(3);
%  r12 = qmat_flattop(5);
%  r21 = qmat_flattop(7);
%  r22 = qmat_flattop(9);
  %
  R = [r11 r12; r21 r22];
  B0_evol(:,:,n+1) = R*B0_evol(:,:,n)*R';
end% for

% calculate max deviation (beta beat) 
b_mat;
b_max_flattop = max(B0_evol(1,1,ceil(end/2):end));
b_min_flattop = min(B0_evol(1,1,ceil(end/2):end));
b_beat = max(abs(b_mat - b_max_flattop), abs(b_mat - b_min_flattop));
b_beat_rel = b_beat / b_mat;


% extend z and np to flat top
z = [z (z+max(z))];
np = [np max(np)*ones(1, length(np))];



if(do_plot)
beta_evol = squeeze(B0_evol(1,1,:));
  set(0,'defaultaxesfontsize',24);
[AX,H1,H2] = plotyy(z, np/1e6, z, beta_evol(1:end-1));

% temp for Yuri
%NN= 66;
%s = z(1:NN)'
%n_p = np(1:NN)'/1e6
%save -ascii ramp.asc s n_p
%hold on;
%plot(z, beta_evol(1:end-1));
grid on;
xlabel('s [m]');
set(get(AX(1),'Ylabel'),'String','n_p [cm^{-3}]'); 
set(get(AX(2),'Ylabel'),'String','\beta [m]');
set(get(AX(1),'Ylabel'),'FontSize',20);
set(get(AX(2),'Ylabel'),'FontSize',20);
end% if
