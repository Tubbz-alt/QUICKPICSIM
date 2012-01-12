% convert elegant particle matrix to quickpic particle matrix
function pp_new=my_qp2ele(pp)
SI_c = 299792458;
SI_em = 9.1093897e-31;
SI_e = 1.60217733e-19;
pp_new(:,1) = pp(:,1)/1e6; % [m]
pp_new(:,2) = pp(:,4)/1e6; % [rad]
pp_new(:,3) = pp(:,2)/1e6; % [m]
pp_new(:,4) = pp(:,5)/1e6; % [rad]
pp_new(:,5) = pp(:,3)/1e6 / SI_c; % [s]
pp_new(:,6) = pp(:,6)*1e9 / (SI_em*SI_c^2/SI_e); % [me * c], or, the Lorentz factor
%  use the following lines for conversion to elegant (EA) instead
%pp_new(:,5) = pp(:,3)/1e6; % [m]
%pp_new(:,6) = pp(:,6)*1e9; % [eV]
