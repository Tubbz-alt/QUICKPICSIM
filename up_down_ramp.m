
% "less /tmp/lens_np.txt" after running calc_ramp.m twice
%
n_down = '1.000e+00,1.000e+00,1.000e+00,9.980e-01,9.956e-01,9.926e-01,9.890e-01,9.845e-01,9.789e-01,9.722e-01,9.640e-01,9.540e-01,9.420e-01,9.275e-01,9.102e-01,8.896e-01,8.654e-01,8.373e-01,8.048e-01,7.680e-01,7.268e-01,6.815e-01,6.327e-01,5.811e-01,5.278e-01,4.738e-01,4.205e-01,3.689e-01,3.201e-01,2.748e-01,2.336e-01,1.968e-01,1.643e-01,1.362e-01,1.120e-01,9.141e-02,7.409e-02,5.960e-02,4.756e-02,3.760e-02,2.939e-02,2.266e-02,1.715e-02,1.265e-02,8.977e-03,5.993e-03,3.567e-03,1.598e-03,1.000e-10,1.000e-10,1.000e-10';
s_down = '0.000e+00,3.000e+05,3.021e+05,3.043e+05,3.064e+05,3.085e+05,3.106e+05,3.128e+05,3.149e+05,3.170e+05,3.191e+05,3.213e+05,3.234e+05,3.255e+05,3.277e+05,3.298e+05,3.319e+05,3.340e+05,3.362e+05,3.383e+05,3.404e+05,3.426e+05,3.447e+05,3.468e+05,3.489e+05,3.511e+05,3.532e+05,3.553e+05,3.574e+05,3.596e+05,3.617e+05,3.638e+05,3.660e+05,3.681e+05,3.702e+05,3.723e+05,3.745e+05,3.766e+05,3.787e+05,3.809e+05,3.830e+05,3.851e+05,3.872e+05,3.894e+05,3.915e+05,3.936e+05,3.957e+05,3.979e+05,4.000e+05,4.021e+05,4.521e+05';

n_up = '1.000e-10,1.000e-10,1.598e-03,3.567e-03,5.993e-03,8.977e-03,1.265e-02,1.715e-02,2.266e-02,2.939e-02,3.760e-02,4.756e-02,5.960e-02,7.409e-02,9.141e-02,1.120e-01,1.362e-01,1.643e-01,1.968e-01,2.336e-01,2.748e-01,3.201e-01,3.689e-01,4.205e-01,4.738e-01,5.278e-01,5.811e-01,6.327e-01,6.815e-01,7.268e-01,7.680e-01,8.048e-01,8.373e-01,8.654e-01,8.896e-01,9.102e-01,9.275e-01,9.420e-01,9.540e-01,9.640e-01,9.722e-01,9.789e-01,9.845e-01,9.890e-01,9.926e-01,9.956e-01,9.980e-01,1.000e+00,1.000e+00,1.000e+00';
s_up = '0.000e+00,2.128e+03,4.255e+03,6.383e+03,8.511e+03,1.064e+04,1.277e+04,1.489e+04,1.702e+04,1.915e+04,2.128e+04,2.340e+04,2.553e+04,2.766e+04,2.979e+04,3.191e+04,3.404e+04,3.617e+04,3.830e+04,4.043e+04,4.255e+04,4.468e+04,4.681e+04,4.894e+04,5.106e+04,5.319e+04,5.532e+04,5.745e+04,5.957e+04,6.170e+04,6.383e+04,6.596e+04,6.809e+04,7.021e+04,7.234e+04,7.447e+04,7.660e+04,7.872e+04,8.085e+04,8.298e+04,8.511e+04,8.723e+04,8.936e+04,9.149e+04,9.362e+04,9.574e+04,9.787e+04,1.000e+05,1.021e+05,4.021e+05';

%q
%
% put together to one profile

% text, for quickpic
n_val = 10;
n_tot = [n_up(1:end-n_val*2+1)  n_down(1+n_val*2:end-n_val*2)];
s_tot = [s_up(1:end-n_val*2+1)  s_down(1+n_val*2:end-n_val*2)];

% auto convert to num, for matlab processing
n = eval(n_tot) % todo ... (for now manually...)

s_ramp = [ 0        2128        4255        6383        8511       10640       12770       14890       17020       19150       21280       23400       25530       27660       29790       31910       34040       36170       38300       40430       42550       44680       46810       48940       51060       53190  55320       57450       59570       61700       63830       65960       68090       70210       72340       74470       76600       78720       80850       82980       85110       87230       89360       91490       93620       95740       97870      100000      302100      304300      306400      308500   310600      312800      314900      317000      319100      321300      323400      325500      327700      329800      331900      334000      336200      338300      340400      342600      344700      346800      348900      351100      353200      355300      357400      359600      361700      363800   366000      368100      370200      372300      374500      376600      378700      380900      383000      385100      387200      389400      391500      393600      395700      397900      400000];
n_ramp = [0.0000    0.0000    0.0016    0.0036    0.0060    0.0090    0.0126    0.0171    0.0227    0.0294    0.0376    0.0476    0.0596    0.0741    0.0914    0.1120    0.1362    0.1643    0.1968    0.2336    0.2748    0.3201    0.3689    0.4205    0.4738    0.5278    0.5811    0.6327    0.6815    0.7268    0.7680    0.8048    0.8373    0.8654    0.8896    0.9102    0.9275    0.9420    0.9540    0.9640    0.9722    0.9789    0.9845    0.9890    0.9926    0.9956    0.9980    1.0000    1.0000    0.9980    0.9956    0.9926    0.9890    0.9845    0.9789    0.9722    0.9640    0.9540    0.9420    0.9275    0.9102    0.8896    0.8654    0.8373    0.8048    0.7680    0.7268    0.6815    0.6327    0.5811    0.5278    0.4738    0.4205    0.3689    0.3201    0.2748    0.2336    0.1968    0.1643    0.1362    0.1120    0.0914    0.0741    0.0596    0.0476    0.0376    0.0294    0.0227    0.0171    0.0126    0.0090    0.0060    0.0036    0.0016    0.0000];

s_ft = [0 150000 300000];
n_ft = [1 1 1];

% down ramp at 70 :
s_ramp(end/2:end) = s_ramp(end/2:end)+4e5;
s_ft(end/2:end) = s_ft(end/2:end)+4e5;


% plugging back for QuickPIC (manually...)
 Density_Variation=.true.
 Density_Variation_NSec=95
 Density_Variation_Fs(1:95) =1.000e-10,1.000e-10,1.598e-03,3.567e-03,5.993e-03,8.977e-03,1.265e-02,1.715e-02,2.266e-02,2.939e-02,3.760e-02,4.756e-02,5.960e-02,7.409e-02,9.141e-02,1.120e-01,1.362e-01,1.643e-01,1.968e-01,2.336e-01,2.748e-01,3.201e-01,3.689e-01,4.205e-01,4.738e-01,5.278e-01,5.811e-01,6.327e-01,6.815e-01,7.268e-01,7.680e-01,8.048e-01,8.373e-01,8.654e-01,8.896e-01,9.102e-01,9.275e-01,9.420e-01,9.540e-01,9.640e-01,9.722e-01,9.789e-01,9.845e-01,9.890e-01,9.926e-01,9.956e-01,9.980e-01,1.000e+00,1.000e+00,9.980e-01,9.956e-01,9.926e-01,9.890e-01,9.845e-01,9.789e-01,9.722e-01,9.640e-01,9.540e-01,9.420e-01,9.275e-01,9.102e-01,8.896e-01,8.654e-01,8.373e-01,8.048e-01,7.680e-01,7.268e-01,6.815e-01,6.327e-01,5.811e-01,5.278e-01,4.738e-01,4.205e-01,3.689e-01,3.201e-01,2.748e-01,2.336e-01,1.968e-01,1.643e-01,1.362e-01,1.120e-01,9.141e-02,7.409e-02,5.960e-02,4.756e-02,3.760e-02,2.939e-02,2.266e-02,1.715e-02,1.265e-02,8.977e-03,5.993e-03,3.567e-03,1.598e-03,1.000e-10
 Density_Variation_s(1:95) = 0.000e+00,2.128e+03,4.255e+03,6.383e+03,8.511e+03,1.064e+04,1.277e+04,1.489e+04,1.702e+04,1.915e+04,2.128e+04,2.340e+04,2.553e+04,2.766e+04,2.979e+04,3.191e+04,3.404e+04,3.617e+04,3.830e+04,4.043e+04,4.255e+04,4.468e+04,4.681e+04,4.894e+04,5.106e+04,5.319e+04,5.532e+04,5.745e+04,5.957e+04,6.170e+04,6.383e+04,6.596e+04,6.809e+04,7.021e+04,7.234e+04,7.447e+04,7.660e+04,7.872e+04,8.085e+04,8.298e+04,8.511e+04,8.723e+04,8.936e+04,9.149e+04,9.362e+04,9.574e+04,9.787e+04,1.000e+05,3.021e+05,3.043e+05,3.064e+05,3.085e+05,3.106e+05,3.128e+05,3.149e+05,3.170e+05,3.191e+05,3.213e+05,3.234e+05,3.255e+05,3.277e+05,3.298e+05,3.319e+05,3.340e+05,3.362e+05,3.383e+05,3.404e+05,3.426e+05,3.447e+05,3.468e+05,3.489e+05,3.511e+05,3.532e+05,3.553e+05,3.574e+05,3.596e+05,3.617e+05,3.638e+05,3.660e+05,3.681e+05,3.702e+05,3.723e+05,3.745e+05,3.766e+05,3.787e+05,3.809e+05,3.830e+05,3.851e+05,3.872e+05,3.894e+05,3.915e+05,3.936e+05,3.957e+05,3.979e+05,4.000e+05


 