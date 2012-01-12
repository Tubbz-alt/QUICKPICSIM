%
% write quickpic data to elegant ascii sdds format
% NB: charge is not written/transferred!
%
% E. Adli, Jan 11, 2012
%
function my_write_qp2elefile(pp_qp, filename);
SI_c = 299792458;
pp = my_qp2ele(pp_qp);
addpath('/Users/eadli/Dropbox/SLAC/litrack/li2elegant');
gam0 = mean(pp(:,6));
% diag matrix is convert to req input for "write_elegant_file()"
% (different from ele format)
pp_2write = pp;
pp_2write(:,5) = pp(:,5) * SI_c;
pp_2write(:,6) = pp(:,6) / gam0 -1;
write_elegant_file(pp_2write, gam0, filename, 'e', 6);
