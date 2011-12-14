% input: matrix, row 1: E [GeV],  row 2, 5 x, xp,  row 3,6 y, yp [m]
function E=my_calc_emittance(B)
   Cxx=sqrt(det(cov(B(:,[2,5]))));
    Cyy=sqrt(det(cov(B(:,[3,6]))));
    E=[Cxx,Cyy]*mean(B(:,1))/0.000510998903076601;
    