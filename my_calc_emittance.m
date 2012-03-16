% input: matrix, row 1: E [GeV],  row 2, 5 x, xp,  row 3,6 y, yp [m]
function E=my_calc_emittance(B)

%size(B,2);
  
if(size(B,2) <= 7)
   Cxx=sqrt(det(cov(B(:,[2,5]))));
    Cyy=sqrt(det(cov(B(:,[3,6]))));
    E=[Cxx,Cyy]*mean(B(:,1))/0.000510998903076601;
    % placet units to SI
    E = E /1e5;
else
    B(:,8:13)=0; % keep in case we want to calc for macro particle clouds as well
    sigma_w = sum(B(:,7));
    if( sigma_w > 0.0 )
    mu_x = B(:,7)'*B(:,2) / sigma_w;
    mu_xp = B(:,7)'*B(:,5) / sigma_w;
    mu_y = B(:,7)'*B(:,3) / sigma_w;
    mu_yp = B(:,7)'*B(:,6) / sigma_w;
    sigma_yy = B(:,7)'*(B(:,11)+((B(:,3)-mu_y).*(B(:,3)-mu_y)));
    sigma_yyp = B(:,7)'*(B(:,12)+((B(:,3)-mu_y).*(B(:,6)-mu_yp)));
    sigma_ypyp = B(:,7)'*(B(:,13)+((B(:,6)-mu_yp).*(B(:,6)-mu_yp)));
    sigma_xx = B(:,7)'*(B(:,8)+((B(:,2)-mu_x).*(B(:,2)-mu_x)));
    sigma_xxp = B(:,7)'*(B(:,9)+((B(:,2)-mu_x).*(B(:,5)-mu_xp)));
    sigma_xpxp = B(:,7)'*(B(:,10)+((B(:,5)-mu_xp).*(B(:,5)-mu_xp)));
      S(1:2,1:2) = [sigma_xx  sigma_xxp; sigma_xxp sigma_xpxp ]/ sigma_w;
      S(3:4,3:4) = [sigma_yy  sigma_yyp; sigma_yyp sigma_ypyp ]/ sigma_w;
      mu_E = B(:,7)'*B(:,1) / sigma_w;
      E = [sqrt(det(S(1:2, 1:2)))  sqrt(det(S(3:4, 3:4))) ] *(mu_E/0.000510998903076601);
    else
      E = [0 0];
    end% if
    % placet units to SI
    E = E / 1e5;
end% if
