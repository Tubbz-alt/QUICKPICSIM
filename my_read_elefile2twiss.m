
% read twiss data from elegant sdds format (bin or ascii)
%
% E. Adli, Jan 11, 2012
%
function [twiss, param] = my_read_elefile2twiss(filename);

elexedir = '/Users/eadli/Dropbox/SLAC/elegant/bin/';
syscmd= [elexedir 'sdds2stream ' filename ' -columns=s,betax,alphax,betay,alphay,etax,etay' ' > /tmp/twiss.txt'];
system(syscmd);
syscmd= [elexedir 'sdds2stream ' filename ' -parameters=nux,nuy,dnux/dp,dnuy/dp' ' > /tmp/param.txt'];
system(syscmd);
twiss = load('/tmp/twiss.txt');
param = load('/tmp/param.txt');
