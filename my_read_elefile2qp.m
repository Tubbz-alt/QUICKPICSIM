%
% read quickpic data from elegant sdds format (bin or ascii)
% NB: charge is not written/transferred!
%
% E. Adli, Jan 11, 2012
%
function pp_qp = my_read_elefile2qp(filename);

elexedir = '/Users/eadli/Dropbox/SLAC/elegant/bin/';
syscmd= [elexedir 'sdds2stream ' filename ' -columns=x,xp,y,yp,t,p' ' > /tmp/data.txt'];
system(syscmd);
pp = load('/tmp/data.txt');
%pp(:,5) =  (pp(:,5)-mean(pp(:,5)));
pp_qp = my_ele2qp(pp);
