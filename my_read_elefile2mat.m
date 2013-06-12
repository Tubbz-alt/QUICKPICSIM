
% read beamline end-to-end matrix from elegant sdds format (bin or ascii)
%
% E. Adli, Jan 11, 2012
%
function [M_cum, M] = my_read_elefile2mat(filename);

elexedir = '/Users/eadli/Dropbox/SLAC/elegant/bin/';
syscmd= [elexedir 'sdds2stream ' filename ' -columns=R11,R12,R13,R14,R15,R16,R21,R22,R23,R24,R25,R26,R31,R32,R33,R34,R35,R36,R41,R42,R43,R44,R45,R46,R51,R52,R53,R54,R55,R56,R61,R62,R63,R64,R65,R66' ' > /tmp/data3.txt'];
system(syscmd);
el_all = load('/tmp/data3.txt');
M_cum = eye(6,6);
%

% the following stores all CUMULATIVE submatrices... (elegant
% stores cum matrices from 1 to n, not from n-1 to n.)
%   (I kept the code below in case I need to get sub-matrices later)
%for n=size(el_all,1):size(el_all,1);   % USE THIS TO SPEED UP IF ONLY M_cum is needed
for n=1:size(el_all,1);
%  el_all(n, 0*6+[1:6])
  M(1,:,n) = el_all(n, 0*6+[1:6]);
  M(2,:,n) = el_all(n, 1*6+[1:6]);
  M(3,:,n) = el_all(n, 2*6+[1:6]);
  M(4,:,n) = el_all(n, 3*6+[1:6]);
  M(5,:,n) = el_all(n, 4*6+[1:6]);
  M(6,:,n) = el_all(n, 5*6+[1:6]);
  %M_cum = M(:,:,n) * M_cum;
end% for
M_cum = M(:,:,n);

%
% get total direcly from elegant
%
%syscmd= [elexedir 'sdds2stream ' '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.fin' ' -parameter=R11,R12,R13,R14,R15,R16,R21,R22,R23,R24,R25,R26,R31,R32,R33,R34,R35,R36,R41,R42,R43,R44,R45,R46,R51,R52,R53,R54,R55,R56' ' > /tmp/data.txt'];
%system(syscmd);
