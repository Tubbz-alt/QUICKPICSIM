
%
% osiris 
%
%data = load('~/temp/test.dat');
%data = load('~/temp/beforeramp.dat'); % Oct 2012
x = data(1:7:end);
y = data(2:7:end);
z = data(3:7:end);
xp = data(4:7:end);
yp = data(5:7:end);
E = data(6:7:end);
q = data(7:7:end);
% fixing offsets
%x = x - mean(x);
%y = y - mean(y);
%xp = xp - mean(x);
%yp = yp - mean(y);
pp_full = [x y z xp yp E];
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam([pp_full q], [1 1 1 0 1 0]);

% downsampling
dist_red = my_reduce_dist([pp_full q], 1, 1504);
%dist_red = [pp_full q];
pp_os = dist_red(:,1:6);
q_red = dist_red(:,7);
filename = '/Users/eadli/Dropbox/SLAC/quickpic/myData/dist_ele.asc';
%my_write_qp2elefile(pp_os, filename);
q = q_red;

% take out low charge tails and create weighted dist
if(1)
[q_sort, I] = sort(q); % sort according to q
cutfrac = 5/6;
charge_pers = sum(q_sort(1:round(end*cutfrac)) / sum(q)) * 100
n_new_particles = sum(sort(q_sort(end*cutfrac:end)) / min(q_sort(end*cutfrac:end)));
q_red_scaled = (sort(q_sort(round(end*cutfrac:end))) / min(q_sort(round(end*cutfrac:end))));
%semilogy(q_red_scaled)
pp_red = pp_os(I, :);
pp_red = pp_red(round(end*cutfrac):end, :);
%
[sigx, sigy, sigxp, sigyp, gauss_sigx, gauss_sigy, gauss_sigxp, gauss_sigyp, emnx, emny, betax, betay, alphax, alphay, muE, sigEE, corzp, corzx, corzy] = my_ana_beam(pp_red, [1 1 1 0 1 0]);
% multiply number of particle for real scaling)
end% if

%pp_mult = my_mult_part(pp_red, round(q_red_scaled));
