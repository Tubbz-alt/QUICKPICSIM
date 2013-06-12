%
% QuickPIC Matlab rpinput generation example script
% E. Adli, Dec 13, 2011
%
% Last update: E. Adli, Jan 15, 2013
%

clear all;
clf;

% update working dir to your QUICKPICSIM folder
working_dir = '/Users/eadli/Dropbox/SLAC/quickpic/QUICKPICSIM/'; eval(['run ' working_dir 'my_SI_params.m']); % import my standard SI constants

rpinput_template_file = [working_dir 'rpinput/rpinput_template'];
rpinput_output_file =  [working_dir 'rpinput/rpinput'];

%
% INPUT TO RPINPUT
%



%
% hosing scen (oct 2012)
%
% from logbook, sx=35, sy=31.  Est: sz=30um
%

%
if(0)
%
ex = 600.0e-6;
ey = 200.0e-6;
E0 = 20;
n0 = 2.66e17;
%
plasma_s_prop = 0.28; % [m]  - propagation length if the beam into the plasma
plasma_density = n0; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.1278;
%lens_demag = 1;
%sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in
                                     % beta and alpha                                    
sigma_x = 35e-6; 
sigma_y = 31e-6;
sigma_z = 30e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = ex;
emitt_y = ey
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if



%
% equilibrum of beam sizes  (Oct 2012)
%
% 2) IDEAL LIN REG EXAMPLE, with flat top (and pre-ion) - with nb/n0 ~ 1/100
%

%
if(0)
%
ex = 600.0e-6;
ey = 200.0e-6;
E0 = 20;
n0 = 2.5e17;
%
plasma_s_prop = 0.30; % [m]  - propagation length if the beam into the plasma
plasma_density = n0; % /cm^3
plasma_PREION=1; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.1278;
%lens_demag = 1;
%sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in
                                     % beta and alpha                                    
sigma_x = 80e-6; % use this for directly into flat top
sigma_y = sigma_x;
sigma_z = 80e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = ex;
emitt_y = ey
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if



%
% equilibrum of beam sizes  (Oct 2012)
%
% 2) REALISTIC E200 2012 example, with flat top (and field ion)
%

%
if(0)
%
ex = 600.0e-6;
ey = 200.0e-6;
E0 = 20;
n0 = 2.5e17;
%
plasma_s_prop = 0.30; % [m]  - propagation length if the beam into the plasma
plasma_density = n0; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.1278;
%lens_demag = 1;
%sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in
                                     % beta and alpha                                    
sigma_x = 40e-6; % use this for directly into flat top
sigma_y = sigma_x;
sigma_z = 40e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = ex;
emitt_y = ey
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if




%
% equilibrum of beam sizes  (Oct 2012)
%
% 1) REALISTIC EXAMPLE WITH RAMP (need to be manually added)
%

%
if(0)
%
ex = 1200.0e-6;
ey = 300.0e-6;
E0 = 20;
n0 = 2.5e17;
%
plasma_s_prop = 0.45; % [m]  - propagation length if the beam into the plasma
plasma_density = n0; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.1278;
%lens_demag = 1;
%sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in
                                     % beta and alpha                                    
sigma_x = 20e-6; % use this for directly into flat top
sigma_y = sigma_x;
sigma_z = 20e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = ex;
emitt_y = ey
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if






%
% linear collider  (Sep-Oct 2012)
%

%  1e17 / cm3
%
if(0)
%
ex = 2e-6;
ey = 0.050e-6;
E0 = 25;
n0 = 1e17;
%
plasma_s_prop = 1.0; % [m]  - propagation length if the beam into the plasma
plasma_density = n0; % /cm^3
plasma_PREION=1; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 1; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.1278;
%lens_demag = 1;
%sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in
                                     % beta and alpha                                    
sigma_x = 20e-6; % use this for directly into flat top
sigma_y = sigma_x;
sigma_z = 20e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = ex;
emitt_y = ey
charge = 1e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if





%
% plasma ramp project - take two  (Sep-Oct 2012)
%

%  1e17 / cm3
%
if(0)
%
bx_mat = 1e-2;
ex = 100e-6;
E0 = 20;
sz_mat = sqrt(ex/E0*.511e-3 * bx_mat)
kp = sqrt(2*(E0/.511e-3)) / bx_mat;
n0 = (kp*SI_c/SI_e)^2 * SI_em * SI_eps0 / 1e6
%
plasma_s_prop = 0.45; % [m]  - propagation length if the beam into the plasma
plasma_density = 2.211e16; % /cm^3
plasma_PREION=1; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
lens_demag = 0.1278;
%lens_demag = 1;
sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in
                                     % beta and alpha                                    
%sigma_x = sz_mat; % use this for directly into flat top
sigma_y = sigma_x;
sigma_z = 20e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = ex;
emitt_y = emitt_x;
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if






%
%  ionization matrix
%
if(0)
% S
plasma_s_prop = 0.4; % [m]  - plasma propagation length
plasma_density = 2.5e17; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=0; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % override sigma_x, sigma_y with matched counterparts
emitt_match =0; % override emitt_x, emitt_y with matched counterparts
mean_E = 39139*0.511*1e6;
sigma_E_E = 0;
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
sigma_x = 20e-6;
sigma_y = sigma_x;
sigma_z = 40e-6;
% this will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = [300.0e-6];
emitt_y = [30.0e-6];
%n_b = 4 * plasma_density;
%charge = n_b * 1e6* sqrt(2*pi*sigma_x^2)*sqrt(2*pi*sigma_y^2)*sqrt(2*pi*sigma_z^2) * SI_e;
charge = 3e-9;
tilt_angle = 0.0;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';
my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match);
end% if



%
% SAREC S-M proposal - gaussian
%
if(0)
% S
plasma_s_prop = 0.1; % [m]  - plasma propagation length
plasma_density = 2.3e17; % /cm^3
plasma_PREION=1; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % override sigma_x, sigma_y with matched counterparts
emitt_match =0; % override emitt_x, emitt_y with matched counterparts
mean_E = 39139*0.511*1e6;
sigma_E_E = 0;
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
sigma_x = 10e-6;
sigma_y = sigma_x;
sigma_z = 500e-6;
% this will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = [70.0e-6];
emitt_y = [70.0e-6];
%n_b = 4 * plasma_density;
%charge = n_b * 1e6* sqrt(2*pi*sigma_x^2)*sqrt(2*pi*sigma_y^2)*sqrt(2*pi*sigma_z^2) * SI_e;
charge = 1.6e-9;
tilt_angle = 0.0;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';
my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match);
end% if






%
% realistic 30^3 + ramp; does ramp help us?
%
%
%  1e17 / cm3
%
if(0)
% S
plasma_s_prop = 0.5; % [m]  - propagation length if the beam into the plasma
plasma_density = 1.0e17; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = 39824*0.511*1e6; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.0879; % example
%lens_demag = 0.105217593817975;
lens_demag = 1; % here we want to test the ramp with QP
sigma_x = 30e-6*sqrt(lens_demag);
sigma_y = sigma_x;
sigma_z = 30e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = 708.0*1e-6 / 2; % realist. emit. in x
emitt_y = 708.0*1e-6 / 8; % real. emit. in y
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if





%
% plasma ramp project
%
%
%  1e17 / cm3
%
if(0)
% S
plasma_s_prop = 0.5; % [m]  - propagation length if the beam into the plasma
plasma_density = 1.0e17; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = 39139*0.511*1e6; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
lens_demag = 0.0879; % example
lens_demag = 0.105217593817975;
%lens_demag = 1;
sigma_x = 20e-6*sqrt(lens_demag);
sigma_y = sigma_x;
sigma_z = 30e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = [350.36*1e-6];
emitt_y = [350.36*1e-6*1e-1];
emitt_y = [350.36*1e-6*1e-0];
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if



%
%  3.6e17 / cm3 - full example (flat top + ramp)
%
if(0)
% S
plasma_s_prop = 0.50; % [m]  - propagation length if the beam into the plasma
plasma_density = 3.6e17; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =1; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = 39139*0.511*1e6; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
lens_demag = 0.09
%lens_demag = 1;
sigma_x = 20e-6*sqrt(lens_demag);
sigma_y = sigma_x;
sigma_z = 30e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = [568.61*1e-6*1e-0];
emitt_y = [568.61*1e-6*1e-0];
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.05;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if






%
%  example input - one beam gaussian
%
if(0)
% S
plasma_s_prop = 0.4; % [m]  - propagation length if the beam into the plasma
plasma_density = 1.0e17; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = 39139*0.511*1e6; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
lens_demag = 0.0879; % example
lens_demag = 0.105217593817975;
%lens_demag = 1;
sigma_x = 20e-6*sqrt(lens_demag);
sigma_y = sigma_x;
sigma_z = 30e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = [350.36*1e-6];
emitt_y = [350.36*1e-6*1e-1];
emitt_y = [350.36*1e-6*1e-0];
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if




%
%  example input - two beam gaussian
%
if(0)
% two-beam test
plasma_s_prop = 0.4; % [m]  - plasma propagation length
plasma_density = 1e17; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step)
beam_match = [1 1]; % override sigma_x, sigma_y with matched counterparts
emitt_match = [0 0]; % override emitt_x, emitt_y with matched counterparts
charge = [3e10*SI_e  0.5*3e10*SI_e]; % [C]
mean_E = [48924*0.511*1e6   0.5*48924*0.511*1e6];
sigma_E_E = [0  0];
sigma_x = [3.28e-6   0.5*3.28e-6];
sigma_y = [3.28e-6   0.5*3.28e-6];
sigma_z = [30.0e-6   30.0e-6];
emitt_x = [100.0e-6  100.0e-6];
emitt_y = [100.0e-6  100.0e-6];
beam_z_pos = [150 200];
tilt_x = [0  0; 0  0; 0  0];
tilt_y = [0  0; 0  0; 0  0];
my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, beam_z_pos, tilt_x, tilt_y);
end% if



%
% PWFALC Fermilab init setup
%
if(0)
% two-beam test
plasma_s_prop = 3.0; % [m]  - plasma propagation length
plasma_density = 1e16; % /cm^3
n0 = plasma_density;
omega_p = sqrt(n0*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
plasma_PREION=1; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=0; % 0 : calc wake only (single time-step)
beam_match = [1 1]; % override sigma_x, sigma_y with matched counterparts
emitt_match = [0 0]; % override emitt_x, emitt_y with matched counterparts
charge = [7.5e10*SI_e  0.5*5e10*SI_e]; % [C]
mean_E = [80e9 80e9];
sigma_E_E = [0  0];
sigma_x = [3.28e-6   0.5*3.28e-6];
sigma_y = [3.28e-6   0.5*3.28e-6];
sigma_z = [1/k_p   1/k_p];
emitt_x = [100.0e-6  100.0e-6];
emitt_y = [100.0e-6  100.0e-6];
beam_z_pos = [180 250];
tilt_x = [0  0; 0  0; 0  0];
tilt_y = [0  0; 0  0; 0  0];
my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, beam_z_pos, tilt_x, tilt_y);
end% if






%
% NJP like example
%
if(0)
% S
plasma_s_prop = 0.4; % [m]  - plasma propagation length
plasma_density = 2.3e17; % /cm^3
plasma_PREION=1; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % override sigma_x, sigma_y with matched counterparts
emitt_match =0; % override emitt_x, emitt_y with matched counterparts
mean_E = 39139*0.511*1e6;
sigma_E_E = 0;
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
sigma_x = 10e-6;
sigma_y = sigma_x;
sigma_z = 500e-6;
% this will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = [70.0e-6];
emitt_y = [70.0e-6];
%n_b = 4 * plasma_density;
%charge = n_b * 1e6* sqrt(2*pi*sigma_x^2)*sqrt(2*pi*sigma_y^2)*sqrt(2*pi*sigma_z^2) * SI_e;
charge = 1.6e-9;
tilt_angle = 0.01;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

% creat gauss
n_sigma = 3.5;
N = 22;
histZ = linspace(0, sigma_z*n_sigma*2, N) - sigma_z*n_sigma;
histZcount = gauss(histZ, 0, sigma_z*1e0);
% cut half for gauss
%histZcount(1:end/2) = 0;

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y, histZ, histZcount);
end% if






%
%  FACET two-bunch init set-up
%
if(0)
% two-beam test
plasma_s_prop = 1.0; % [m]  - plasma propagation length
plasma_density = 3e16; % /cm^3
n0 = plasma_density;
omega_p = sqrt(n0*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
plasma_PREION=0; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step)
beam_match = [0 0]; % override sigma_x, sigma_y with matched counterparts
emitt_match = [0 0]; % override emitt_x, emitt_y with matched counterparts
charge = [1.5e10*SI_e  0.5e10*SI_e]; % [C]
mean_E = [20e9 20e9];
sigma_E_E = [0  0];
sigma_x = [20.0e-6   20.0e-6];
sigma_y = [20.0e-6   20.0e-6];
sigma_z = [30e-6   20e-6];
emitt_x = [300.0e-6  300.0e-6];
emitt_y = [30.0e-6  30.0e-6];
beam_z_pos = [180 330];
tilt_x = [0  0; 0  0; 0  0];
tilt_y = [0  0; 0  0; 0  0];
my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, beam_z_pos, tilt_x, tilt_y);
end% if






%
% PWFALC - NEW OPTIM
%
if(0)
% two-beam test
T = 1.5;
plasma_s_prop = 2.2; % [m]  - plasma propagation length
plasma_density = 2e16; % /cm^3
n0 = plasma_density;
omega_p = sqrt(n0*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
plasma_PREION=1; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=0; % 0 : calc wake only (single time-step)
beam_match = [1 1]; % override sigma_x, sigma_y with matched counterparts
emitt_match = [0 0]; % override emitt_x, emitt_y with matched counterparts
charge = [4.5e10*SI_e  1.0e1*SI_e]; % [C]
mean_E = [17e9 25e9];
sigma_E_E = [0  0];
sigma_x = [3.28e-6   0.5*3.28e-6]; % N/a
sigma_y = [3.28e-6   0.5*3.28e-6]; % N/a
sigma_z = [1/k_p   1/k_p];
emitt_x = [2.0e-6  2.0e-6];
emitt_y = [2.0e-6  2.0e-6];
beam_z_pos = [4*(1/k_p)*1e6  (4+6)*(1/k_p)*1e6];
tilt_x = [0  0; 0  0; 0  0];
tilt_y = [0  0; 0  0; 0  0];
my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, beam_z_pos, tilt_x, tilt_y);
end% if









%
% FACET-II  (Jan 2013)
%
if(0)
% two-beam test
plasma_s_prop = 2.2; % [m]  - plasma propagation length
plasma_density = 1e15; % /cm^3
n0 = plasma_density;
omega_p = sqrt(n0*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
plasma_PREION=1; % 0 : non-ionized plasma (Nneutrals=1) 
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=0; % 0 : calc wake only (single time-step)
beam_match = [1 1]; % override sigma_x, sigma_y with matched counterparts
emitt_match = [0 0]; % override emitt_x, emitt_y with matched counterparts
charge = [3.0  0.3]*1e-9; % [C]
mean_E = [10e9 0.1e9];
sigma_E_E = [0  0];
sigma_x = [3.28e-6   0.5*3.28e-6]; % N/a
sigma_y = [3.28e-6   0.5*3.28e-6]; % N/a
sigma_z = [30e-6   30e-6];
emitt_x = [10.0e-6  1.0e-6];
emitt_y = [10.0e-6  1.0e-6];
beam_z_pos = [4*(1/k_p)*1e6  (4+4)*(1/k_p)*1e6];
tilt_x = [0  0; 0  0; 0  0];
tilt_y = [0  0; 0  0; 0  0];
my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, beam_z_pos, tilt_x, tilt_y);
end% if






%
% CTM studies (jan 2013)
%

%
if(1)
%
ex = 600.0e-6; % n/a
ey = 200.0e-6; % n/a
E0 = 20;
n0 = 1.0e17;
%
plasma_s_prop = 0.55; % [m]  - propagation length if the beam into the plasma
plasma_density = n0; % /cm^3
plasma_PREION=1; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =1; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.1278;
%lens_demag = 1;
%sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in beta and alpha)                                    
%sigma_x = 3.5355e-6; 
sigma_x = 5e-6 * sqrt(2); 
%sigma_x = 5e-6 / sqrt(16384); 
%sigma_x = 5e-6 / sqrt(262144); 
sigma_y = sigma_x;
sigma_z = 30e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = ex;
emitt_y = ey
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.05;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if











%
% Ken/UCLA short bunch "proposal"
%
if(1)
%
ex = 100e-6;
ey = 10e-6;
E0 = 20;
n0 = 4e17;
%
plasma_s_prop = 0.5; % [m]  - propagation length if the beam into the plasma
plasma_density = n0; % /cm^3
plasma_PREION=1; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.1278;
%lens_demag = 1;
%sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in
                                     % beta and alpha                                    
sigma_x = 30e-6; % use this for directly into flat top
sigma_y = sigma_x;
sigma_z = 30e-6;
% 
emitt_x = ex;
emitt_y = ey
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.00;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if







%
% My QP run
%
if(1)
%
ex = 71.6e-6;
ey = 7.16e-6;
E0 = 20.35;
n0 = 2.66e17;
%
plasma_s_prop = 0.28; % [m]  - propagation length if the beam into the plasma
plasma_density = n0; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 37; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 0; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = E0*1e9; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p
%lens_demag = 0.1278;
%lens_demag = 1;
%sigma_x = sz_mat / sqrt(lens_demag); % use this for ramped (put in
                                     % beta and alpha                                    
sigma_x = 30e-6; % use this for directly into flat top
sigma_y = sigma_x;
sigma_z = 60e-6;
% 
emitt_x = ex;
emitt_y = ey
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.5;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match, -1, tilt_x, tilt_y);
end% if
