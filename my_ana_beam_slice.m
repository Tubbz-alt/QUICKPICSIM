%
% my_ana_beam_slice( slice, myaxis, myrange, s_timestep );
%
% E. Adli, 2012
%
% sliced beam display
function my_ana_beam_slice( slice, myaxis, myrange, s_timestep );

if( ~exist('myrange' ) )
  % don't show slices with very small fraction of charge (outlier behaviour)
  n_frac_acceptance = 1.0; % req at least this percentage of slice
  (slice.N_z / sum(slice.N_z) * 100);
  n_good = (slice.N_z / sum(slice.N_z) * 100) > n_frac_acceptance;
  myrange = min(find(n_good)) : max(find(n_good));
  %n_z_range = min(find(n_good)) : max(find(n_good)) -1
end% if

% <x>    
%  n_take_out_avg = 5:8;
%  plot(slice.z(myrange), slice.mean_x - mean(slice.mean_x(n_take_out_avg)), '-x');
  subplot(2,3,2);
    %
    % <eps_y>
    %
plot(slice.z(myrange), slice.emny(myrange), '-x');
  %xlabel('z [um]');
  ylabel('\epsilon_{N,y} [um]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      axis([myaxis.z_min myaxis.z_max myaxis.emnx_min myaxis.emnx_max]);
    end; end;
  %  
  plot(slice.z(myrange), slice.mean_x(myrange), '-x');
  xlabel('z [um]');
  ylabel('<x> [um]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      axis([myaxis.z_min myaxis.z_max myaxis.mean_x_min myaxis.mean_x_max]);
    end; end;


  %  alphax    
  %subplot(2,3,2);
  %plot(slice.z(myrange), slice.alphax(myrange), '-x');
  %xlabel('z [um]');
  %ylabel('\alpha_x [-]');
  %grid on;
  %if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
  %    axis([myaxis.z_min myaxis.z_max myaxis.alphax_min myaxis.alphax_max]);
  %  end; end;
  %
  subplot(2,3,5);
  plot(slice.z(myrange), slice.sigma_x(myrange), '-x');
  xlabel('z [um]');
  ylabel('\sigma_x [um]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      axis([myaxis.z_min myaxis.z_max myaxis.sigma_x_min myaxis.sigma_x_max]);
    end; end;
  %  

  subplot(2,3,4);
  plot(slice.z(myrange), slice.mean_E(myrange), '-x');
  xlabel('z [um]');
  ylabel('<E> [GeV]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      axis([myaxis.z_min myaxis.z_max myaxis.mean_E_min myaxis.mean_E_max]);
    end; end;
  %  
  subplot(2,3,3);
  plot(slice.z(myrange), slice.emnx(myrange), '-x');
  %xlabel('z [um]');
  ylabel('\epsilon_{N,x} [um]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      axis([myaxis.z_min myaxis.z_max myaxis.emnx_min myaxis.emnx_max]);
    end; end;
  %  

  %
  subplot(2,3,5);
    plot(slice.z(myrange), slice.mean_xp(myrange), '-x');
  xlabel('z [um]');
  ylabel('<xp> [urad]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      axis([myaxis.z_min myaxis.z_max myaxis.mean_xp_min myaxis.mean_xp_max]);
    end; end;
%
  plot(slice.z(myrange), slice.betay(myrange), '-x');
  xlabel('z [um]');
  ylabel('\beta_y [m]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      axis([myaxis.z_min myaxis.z_max myaxis.betax_min myaxis.betax_max]);
    end; end;
  %

  subplot(2,3,6);
  plot(slice.z(myrange), slice.betax(myrange), '-x');
  xlabel('z [um]');
  ylabel('\beta_x [m]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      axis([myaxis.z_min myaxis.z_max myaxis.betax_min myaxis.betax_max]);
    end; end;
  %
  subplot(2,3,1);
  bar(slice.z(myrange), slice.N_z(myrange) / max(slice.N_z(myrange)));
  %xlabel('z [um]');
  ylabel('\lambda_z [arb. units]');
  grid on;
  if(exist('myaxis')), if(isfield(myaxis, 'z_min')), 
      myaxis_N = axis;
      axis([myaxis.z_min myaxis.z_max myaxis_N(3) myaxis_N(4)]);
    end; end;
  %  
  
  if(exist('s_timestep')),
    title(['s=' num2str(s_timestep*100, 3) ' [cm].' ]);
  end% if