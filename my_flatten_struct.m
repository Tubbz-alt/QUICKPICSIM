function data = my_flatten_struct(qp, myvar)

for n=1:size(qp,2),
  for m=1:size(qp(1).PP,2),
  eval(['data(n, m) = qp(n).PP(m).' myvar ';']);
  end% for
end% for

return;

