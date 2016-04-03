function [H_err, H_coverage] = plotPErecovery ( PErecovery )

  set (0, "defaultlinemarkersize", 3);

  SNR_inj = [ PErecovery.SNR_inj ];
  max_x = 1.1 * max ( SNR_inj(:) );
  %% ====================
  H_err = figure(); clf;
  %% ---------- relerr(A)
  subplot ( 3, 2, 1 ); hold on;
  y = [PErecovery.A_relerr];
  plot ( SNR_inj, y, "+" );
  min_y = -1; max_y = 1;
  indsAbove = find ( y > max_y );
  if ( !isempty(indsAbove) ) plot ( SNR_inj(indsAbove), max_y * ones(size(indsAbove)), "^" ); endif
  indsBelow = find ( y < min_y );
  if ( !isempty(indsBelow) ) plot ( SNR_inj(indsBelow), min_y * ones(size(indsBelow)), "v" ); endif
  ylim ( [min_y, max_y ] );
  xlim ( [ 0, max_x ] );
  ylabel ( "relerr(A)" ); grid on;
  %% ---------- relerr(phi0)
  subplot ( 3, 2, 2 ); hold on;
  y = [PErecovery.phi0_relerr];
  plot ( SNR_inj, y, "+" );
  min_y = 0; max_y = 0.5;
  indsAbove = find ( y > max_y );
  if ( !isempty(indsAbove) ) plot ( SNR_inj(indsAbove), max_y * ones(size(indsAbove)), "^" ); endif
  ylim ( [min_y, max_y ] );
  xlim ( [ 0, max_x ] );
  ylabel ( "err(phi0)/2pi" ); grid on;
  %% ---------- relerr(f0)
  subplot ( 3, 2, 3 ); hold on;
  y = [ PErecovery.f0_relerr];
  plot ( SNR_inj, y, "+" );
  min_y = -0.5; max_y = 0.5;
  indsAbove = find ( y > max_y );
  if ( !isempty(indsAbove) ) plot ( SNR_inj(indsAbove), max_y * ones(size(indsAbove)), "^" ); endif
  indsBelow = find ( y < min_y );
  if ( !isempty(indsBelow) ) plot ( SNR_inj(indsBelow), min_y * ones(size(indsBelow)), "v" ); endif
  ylim ( [min_y, max_y ] );
  xlim ( [ 0, max_x ] );
  ylabel ( "relerr(f0)" ); grid on;
  %% ----------
  subplot ( 3, 2, 4 ); hold on;
  y = [ PErecovery.tau_relerr ];
  plot ( SNR_inj, y, "+" );
  min_y = -1; max_y = 1;
  indsAbove = find ( y > max_y );
  if ( !isempty(indsAbove) ) plot ( SNR_inj(indsAbove), max_y * ones(size(indsAbove)), "^" ); endif
  indsBelow = find ( y < min_y );
  if ( !isempty(indsBelow) ) plot ( SNR_inj(indsBelow), min_y * ones(size(indsBelow)), "v" ); endif
  ylim ( [min_y, max_y ] );
  xlim ( [ 0, max_x ] );
  ylabel ( "relerr(tau)" ); grid on;
  %% ----------
  subplot ( 3, 2, 5 ); hold on;
  plot ( SNR_inj, SNR_inj, "-o;SNR-inj;", "color", "green" );
  plot ( SNR_inj, [ PErecovery.SNR_MP ], "+;SNR-MP;" );
  plot ( SNR_inj, [ PErecovery.SNR_ML ], "x;SNR-ML;" );
  xlim ( [ 0, max_x ] );
  ylim ( [ 0, max([ PErecovery.SNR_MP, PErecovery.SNR_ML, SNR_inj ](:) ) ] );
  grid on;
  legend ( "location", "northwest" );
  xlabel ( "SNR-inj" ); ylabel ("SNR");
  %% ----------
  subplot ( 3, 2, 6 ); hold on;
  log10BSG = log10 ( [ PErecovery.BSG ] );
  plot ( SNR_inj, log10BSG, "o" );
  min_y = -2; max_y = 5;
  indsAbove = find ( log10BSG > max_y );
  if ( !isempty(indsAbove) ) plot ( SNR_inj(indsAbove), max_y * ones(size(indsAbove)), "^" ); endif
  indsBelow = find ( log10BSG < min_y );
  if ( !isempty(indsBelow) ) plot ( SNR_inj(indsBelow), min_y * ones(size(indsBelow)), "v" ); endif
  xlim ( [ 0, max_x ] );
  ylim ( [min_y, max_y ] );
  line ( xlim(), [1,1] );
  grid on;
  xlabel ( "SNR-inj" );
  ylabel ( "log10(BSG)");
  %% ----------

  H_coverage = figure();
  perc = [PErecovery.f0_tau_percentile] ( find ( [PErecovery.BSG] > 1 ) ); %% restrict 'coverage' tests to things we'd actually call 'detections'
  perc = sort ( perc );
  Nn = length(perc);
  cdf = [1:Nn] / Nn;
  hold on;
  stairs ( perc, cdf );
  plot ( [0,1], [0,1], "-", "color", "green" );
  legend ( "measured", "exact" );
  grid on;
  xlabel ( "nominal coverage" );
  ylabel ( "posterior-coverage (BSG>1)" );
  legend ( "location", "northwest" );

  return;

endfunction
