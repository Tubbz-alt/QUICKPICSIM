% convert elegant particle matrix to quickpic particle matrix
function pp_new=my_qp2ele(pp)
working_dir = '/Users/eadli/Dropbox/SLAC/quickpic/myScripts/'; eval(['run ' working_dir 'my_SI_params.m']); % import my standard SI constants
pp_new(:,1) = pp(:,1)/1e6; % [m]
pp_new(:,2) = pp(:,4)/1e6; % [rad]
pp_new(:,3) = pp(:,2)/1e6; % [m]
pp_new(:,4) = pp(:,5)/1e6; % [rad]
pp_new(:,5) = pp(:,3)/1e6 / SI_c; % [s]
pp_new(:,6) = pp(:,6)*1e9 / (SI_em*SI_c^2/SI_e); % [me * c], or, the Lorentz factor
