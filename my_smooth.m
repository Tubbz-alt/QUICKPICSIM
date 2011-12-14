function y = my_smooth(x, n_filter)
   x = reshape(x, 1, length(x));
   n_half = floor((n_filter-1)/2);
   x_pad = [zeros(1,n_half) x zeros(1,n_half)];

   for(n=1:length(x))
     y(n) = mean(x_pad( (n:(n+2*n_half) )));
   end;
%return y;
   
