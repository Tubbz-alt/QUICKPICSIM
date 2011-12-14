%
% QuickPIC Matlab rpinput generation example script
% E. Adli, Dec 13, 2011
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
%  example input - one beam gaussian
%
if(1)
% S
plasma_s_prop = 1.0; % [m]  - propagation length if the beam into the plasma
plasma_density = 5.4e16; % /cm^3
plasma_PREION=0; % 0 : non-ionized plasma 1: pre-ionized plasma
plasma_Z = 3; % atomic number of plasma gas
BEAM_EV=1; % 0 : calc wake only (single time-step):
           % assume: use 8 cores (express) when set to 0, 
           % and 128 cores when set to 1
beam_match = 1; % 1: override sigma_x, sigma_y with matched counterparts, 0: do nothing
emitt_match =0; % 1: override emitt_x, emitt_y with matched counterparts, 0: do nothing
mean_E = 39139*0.511*1e6; % beam mean energy [eV]
sigma_E_E = 0; % 
omega_p = sqrt(plasma_density*1e6* SI_e^2 /SI_em / SI_eps0);  % plasma frequency
lambda_p = SI_c / omega_p * 2*pi; % plasma wavelength
k_p = 2*pi/lambda_p;
sigma_x = 20e-6;
sigma_y = sigma_x;
sigma_z = 20e-6;
% this following will enforce (close to) ideal density
%k_p = sqrt(2) / sigma_z;
%n0 = (k_p*SI_c)^2 * SI_em *SI_eps0 / SI_e^2  / 1e6;
%plasma_density = n0;
% 
emitt_x = [300.0e-6];
emitt_y = [30.0e-6];
charge = 2e10*SI_e; % beam electron charge [C]
tilt_angle = 0.0;
tilt_x = [0 tilt_angle 0]';
tilt_y = [0 0 0]';

my_gen_rpinput(rpinput_template_file, rpinput_output_file, plasma_density, plasma_Z, plasma_PREION, plasma_s_prop, charge, mean_E, sigma_E_E, sigma_x, sigma_y, sigma_z, emitt_x, emitt_y, BEAM_EV, beam_match, emitt_match);
end% if




%
%  example input - two beam gaussian
%
if(0)
% two-beam test
plasma_s_prop = 0.8; % [m]  - plasma propagation length
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
% example input - arbitrary z-profile, one beam
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

