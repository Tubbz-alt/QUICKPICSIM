% convert elegant particle matrix to quickpic particle matrix
function pp_qp=my_ele2qp(pp)
pp_qp(:,1) = pp(:,1)*1e6;
pp_qp(:,2) = pp(:,3)*1e6;
pp_qp(:,3) = pp(:,5)*1e6;
pp_qp(:,4) = pp(:,2)*1e6;
pp_qp(:,5) = pp(:,4)*1e6;
pp_qp(:,6) = pp(:,6)/1e9;
