% convert elegant particle matrix to quickpic particle matrix
function pp_qp=my_ele2qp(pp)
SI_c = 299792458;
SI_em = 9.1093897e-31;
SI_e = 1.60217733e-19;
pp_qp = pp;
pp_qp(:,1) = pp(:,1)*1e6;
pp_qp(:,2) = pp(:,3)*1e6;
pp_qp(:,4) = pp(:,2)*1e6;
pp_qp(:,5) = pp(:,4)*1e6;
pp_qp(:,3) = pp(:,5)*1e6 * SI_c; % [um]
pp_qp(:,6) = pp(:,6)/1e9 * (SI_em*SI_c^2/SI_e); % [GeV]
%  use the following lines for conversion to elegant (EA) instead
%pp_qp(:,3) = pp(:,5)*1e6;
%pp_qp(:,6) = pp(:,6)/1e9;
