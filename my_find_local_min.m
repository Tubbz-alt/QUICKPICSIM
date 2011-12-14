function n_min = my_find_local_min(z_count, n_smooth)

if( nargin < 3 )
  n_smooth = round(length(z_count)/10);
end% if

%z_count_lp = my_SG_filter(z_count, 3, n_smooth);
z_count_lp = my_smooth(z_count, n_smooth);

%plot(z_count, 'b')
%hold on;
%plot(z_count_lp, 'r')
%hold off;

% avoid local mins at edge
n = round(length(z_count) / 10);

% find first peak
while( (z_count_lp(n+1) ) >= z_count_lp(n) )
  n = n + 1;
end% while

% find valley
while( (z_count_lp(n+1) ) <= z_count_lp(n) )
  n = n + 1;
end% while

n_min = n;
