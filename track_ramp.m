%
% get initial beam from quickpic
%

load /Users/eadli/Dropbox/SLAC/quickpic/myData/n3e17_30cm.mat
charge = 3e-9; % NB: charge is not written/transferred!
pp_qp = qp(1).PP(1).BEAM;
pp_qp = my_reduce_dist(pp_qp, 2); % reduce by a factor N for quicker running
%pp_qp(:,6) = 20e0 * (1 + 1e-6*randn(size(pp_qp(:,6)))); % to make monoenergetic beam
%pp_qp(:,6) = pp_qp(:,6) * 100;
filename = '/Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc';
my_write_qp2elefile(pp_qp, filename);




%
% 1) calc and export plasma ramp focusing strength
%

% plasma density for ramp
n0 = 3.6e17;
% ref energy (should be equal to ref energy of rest of simulation)
if( exist('pp_qp') )
  gamma = mean(pp_qp(:,6));
else
  gamma = 20e3/0.511;
end% if
gamma = 20e3/0.511
my_calc_ramp(n0, gamma, 0, 0);

% 
% 2) run track_lens.sh in elesim
%    cd /Users/eadli/Dropbox/SLAC/elegant/FACETSIM/ELESIM
%    sh ./track_lens.sh 
%

% 
% 3) generate and export beam after plasma
%
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet_lensX.out';
pp_X = my_read_elefile2qp(filename);
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet_lensY.out';
pp_Y = my_read_elefile2qp(filename);

pp_XY = [ pp_X(:,1) pp_Y(:,2:3) pp_X(:,4) pp_Y(:,5:6) ];

filename = '/Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc';
my_write_qp2elefile(pp_XY, filename);
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp_XY, [1 1 1 0 1 0 1]);



%
% 4) track beam in dump line
%    cd /Users/eadli/Dropbox/SLAC/elegant/FACETSIM/ELESIM
%    bash
%    source ../SOURCE
%   ./sim_elegant.sh -i /Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc -d ../FACET_OPTICS/facet_v27.4.dynamic.ele   -b FACET -n off   -r ../FACET_OPTICS/R56Params/10.0mmR56.par   -l ../FACET_OPTICS/facet_simple_dump_only.lte    -o ../FACET_WORK/facet.out
%    %%%%./sim_elegant.sh -i /Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc -d ../FACET_OPTICS/facet_v27.4.dynamic.ele   -b FACET -n off   -r ../FACET_OPTICS/R56Params/10.0mmR56.par   -l ../FACET_OPTICS/facet_v27.4.dump_only.lte    -o ../FACET_WORK/facet.out


%
% analyze final beam
%
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.out';
pp_DUMP = my_read_elefile2qp(filename);
%filename_nomag = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.out_nomagnets';
%pp_DUMP = my_read_elefile2qp(filename_nomag);
%[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp_DUMP, [1 1 1 0 1 0 0]);
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy, slice] = my_ana_beam(pp_DUMP, [0 0 0  0 1 0 0]); sigx, sigy
%
% CHECK IMAGING CONDITION
%
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.mat';
[M, M_arr] = my_read_elefile2mat(filename);
m12 = M(1,2) 
m34 = M(3,4) 

%
% for analysis: twiss in plasma lens
%

filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet_lensX.twi';
[twiss_X] = my_read_elefile2twiss(filename);
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet_lensY.twi';
[twiss_Y] = my_read_elefile2twiss(filename);
n_lensend = length(twiss_X);
%ax = sqrt(twiss(end,2) / twiss(n_lensend,2))
%ay = sqrt(twiss(end,4) / twiss(n_lensend,4))
frac_plot = n_lensend / length(twiss_X);
plot(twiss_X(1:end*frac_plot,1), twiss_X(1:end*frac_plot,2), 'b')
hold on;
plot(twiss_X(1:end*frac_plot,1), twiss_Y(1:end*frac_plot,4), 'r')
hold off;
xlabel('s [m]');
ylabel('beta [m]');
legend('x', 'y');
grid on;



%
% Calc ax, ay
%
filename = '/Users/eadli/Dropbox/SLAC/elegant/FACETSIM/FACET_WORK/facet.twi';
[twiss, param] = my_read_elefile2twiss(filename);
%n_lensend = (64*3+2);
n_lensend = 1; %length(twiss);
ax = sqrt(twiss(end,2) / twiss(n_lensend,2))
ay = sqrt(twiss(end,4) / twiss(n_lensend,4))
%plot(twiss(1:end,1), twiss(1:end,2));
frac_plot = n_lensend / length(twiss);
%frac_plot = 6/6;
plot(twiss(1:end*frac_plot,1), twiss(1:end*frac_plot,2), 'b')
hold on;
plot(twiss(1:end*frac_plot,1), twiss(1:end*frac_plot,4), 'r')
hold off;
xlabel('s [m]');
ylabel('beta [m]');
legend('x', 'y');
grid on;


%
%
%
%    p0 = 4 GeV
%E = [1 2
%n_loss = [38032  0 
%    p0 = 8 GeV
%E = [2 3
%n_loss = [19801  0 
%    p0 = 10 GeV
%E = [3 4
%n_loss = [186 0
%    p0 = 20 GeV
%E = [4 5 6
%n_loss = [97120 4038 0
     