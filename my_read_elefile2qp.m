%
% read quickpic data from elegant sdds format (bin or ascii)
% NB: charge is not written/transferred!
%
% E. Adli, Jan 11, 2012
%   first version
% E. Adli, Jun 2, 2013
%   added option for 7th coord (pid)
function [pp_qp, pp_qp7] = my_read_elefile2qp(filename);

elexedir = '/Users/eadli/Dropbox/SLAC/elegant/bin/';
syscmd= [elexedir 'sdds2stream ' filename ' -columns=x,xp,y,yp,t,p,particleID' ' > /tmp/data.txt'];
system(syscmd);
pp = load('/tmp/data.txt');
%pp(:,5) =  (pp(:,5)-mean(pp(:,5)));
if( size(pp,1) > 0)
  pp_qp = my_ele2qp(pp(:,1:6));
  pp_qp7 = my_ele2qp(pp(:,1:7));
else
  pp_qp = [];
  pp_qp7 = [];
end% if