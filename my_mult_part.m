% increase particle dist by multiplying each particle by corresponding number in vector N
function pp_mult=my_mult_part(pp, N)

m_tot = 0;
for(n=1:length(N))
  n_mult = N(n);
  for o=1:n_mult;
    pp_mult(m_tot+o, :) = pp(n, :);
  end% for
  m_tot = m_tot + n_mult;
  n;
end% for
