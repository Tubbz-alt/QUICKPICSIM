% input: matrix, row 2, 5 x, xp,  row 3,6 y, yp
function S=my_calc_twiss(B)
    Cx=cov(B(:,[2,5]));
    Cy=cov(B(:,[3,6]));
    Ex=sqrt(det(Cx));
    Ey=sqrt(det(Cy));
    S=zeros(4,4);
    Cx = Cx / Ex;
    Cy = Cy / Ey;
    S(1:2,1:2)=Cx;
    S(3:4,3:4)=Cy;
