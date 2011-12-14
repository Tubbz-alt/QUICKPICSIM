% Savitzky-Golay optimal filter
%   EA, nov 2010
function y = my_SG_filter(x, SG_order, SG_halfwidth)
  %
  % calc of coefficients
  %
  SG_length = 2*SG_halfwidth + 1;
  for(m=0:SG_order)
    for(j=0:SG_length-1),
       A(j+1, m+1) = (j-SG_halfwidth)^m;
    end% for
  end% for
  for(j=0:SG_length-1), 
    Ae = A(j+1,:);
    SG_out = inv(A'*A)*Ae';
    SG_coeffs(j+1) = SG_out(1);
  end% for
                                                                                                            %
  % filtering
  %
  x = reshape(x, 1, length(x));
  x_pad = [zeros(1,SG_halfwidth) x zeros(1,SG_halfwidth)];

   for(n=1:length(x))
     %y_lp(n) = mean(x_pad( (n:(n+2*SG_halfwidth) )));
     y_SG(n) = sum((x_pad( (n:(n+2*SG_halfwidth) ))).*SG_coeffs);
   end;
   y = y_SG';
  
end% function
  

  
