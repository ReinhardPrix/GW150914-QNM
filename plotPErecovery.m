function [H_err, H_coverage] = plotPErecovery ( PErecovery, SNR_limit = 30 )
  global debugLevel = 1;

  set (0, "defaultlinemarkersize", 3);

  if ( !isempty ( SNR_limit ) )
    useBins = find ( [PErecovery.SNR_inj] < SNR_limit );
    PErecovery = PErecovery(useBins);
  endif

  %% nuke all overflowed BSGs (too loud injections)
  mInf = find ( isinf ( [PErecovery.BSG] ) );
  PErecovery ( mInf ) = [];

  SNR_inj = [ PErecovery.SNR_inj ];

  %% ----- histogram the data over bins for 'data-compression'
  x = SNR_inj;
  x0 = round ( min(x) );
  x1 = round ( max(x) );
  xVals = [ x0 : x1 ];	%% put bin-centers on integer values of SNR-inj
  num_xBins = length(xVals);
  binsL = xVals - 0.5;
  binsR = xVals + 0.5;
  hgrms.A = hgrms.phi0 = hgrms.f0 = hgrms.tau = hgrms.SNR = hgrms.log10BSG = cell ( 1, num_xBins );
  for i = 1:num_xBins
    inds_xBin_i = find ( (x >= binsL(i)) & (x < binsR(i)) );
    if ( isempty ( inds_xBin_i ) ) continue; endif

    hgrms.A{i}.hgrm = Hist ( 1, {"lin", "dbin", 0.01 } );
    hgrms.A{i}.hgrm = addDataToHist ( hgrms.A{i}.hgrm, [PErecovery.A_relerr]'(inds_xBin_i) );

    hgrms.phi0{i}.hgrm = Hist ( 1, {"lin", "dbin", 0.01 } );
    hgrms.phi0{i}.hgrm = addDataToHist ( hgrms.phi0{i}.hgrm, [PErecovery.phi0_relerr]'(inds_xBin_i) );

    hgrms.f0{i}.hgrm = Hist ( 1, {"lin", "dbin", 0.01 } );
    hgrms.f0{i}.hgrm = addDataToHist ( hgrms.f0{i}.hgrm, [PErecovery.f0_relerr]'(inds_xBin_i) );

    hgrms.tau{i}.hgrm = Hist ( 1, {"lin", "dbin", 0.01 } );
    hgrms.tau{i}.hgrm = addDataToHist ( hgrms.tau{i}.hgrm, [PErecovery.tau_relerr]'(inds_xBin_i) );

    hgrms.SNR{i}.hgrm = Hist ( 1, {"lin", "dbin", 0.01 } );
    hgrms.SNR{i}.hgrm = addDataToHist ( hgrms.SNR{i}.hgrm, [PErecovery.SNR_MP]'(inds_xBin_i) );

    hgrms.log10BSG{i}.hgrm = Hist ( 1, {"lin", "dbin", 0.01 } );
    hgrms.log10BSG{i}.hgrm = addDataToHist ( hgrms.log10BSG{i}.hgrm, log10([PErecovery.BSG](inds_xBin_i))' );

  endfor
  max_x = max(binsR);

  %% ====================
  H_err = figure(); clf;
  %% ---------- relerr(A)
  subplot ( 3, 2, 1 ); hold on;
  plotHists ( xVals, hgrms.A );
  min_y = -1; max_y = 1;
  ylim ( [min_y, max_y ] );
  xlim ( [ 0, max_x ] );
  ylabel ( "relerr(A)" );
  legend ( "location", "southeast" );

  %% ---------- relerr(phi0)
  subplot ( 3, 2, 2 ); hold on;
  plotHists ( xVals, hgrms.phi0 );
  min_y = 0; max_y = 0.5;
  ylim ( [min_y, max_y ] );
  xlim ( [ 0, max_x ] );
  ylabel ( "err(phi0)/2pi" );
  legend ( "location", "northeast" );

  %% ---------- relerr(f0)
  subplot ( 3, 2, 3 ); hold on;
  plotHists ( xVals, hgrms.f0 );
  min_y = -0.5; max_y = 0.5;
  ylim ( [min_y, max_y ] );
  xlim ( [ 0, max_x ] );
  ylabel ( "relerr(f0)" );

  %% ---------- relerr(tau)
  subplot ( 3, 2, 4 ); hold on;
  plotHists ( xVals, hgrms.tau );
  min_y = -1; max_y = 1;
  ylim ( [min_y, max_y ] );
  xlim ( [ 0, max_x ] );
  ylabel ( "relerr(tau)" );

  %% ---------- SNR_MP
  subplot ( 3, 2, 5 ); hold on;
  plotHists ( xVals, hgrms.SNR );
  plot ( SNR_inj, SNR_inj, "--;SNR-inj;", "color", "green", "linewidth", 2 );
  xlim ( [ 0, max_x ] );
  %%ylim ( [ 0, 1.1 * max( [PErecovery.SNR_MP] ) ] );
  legend ( "location", "northwest" );
  xlabel ( "SNR-inj" ); ylabel ("SNR-MP");

  %% ---------- log10BSG
  subplot ( 3, 2, 6 ); hold on;
  plotHists ( xVals, hgrms.log10BSG );
  min_y = -2; max_y = 5;
  xlim ( [ 0, max_x ] );
  %%ylim ( [min_y, max_y ] );
  line ( xlim(), [0,0], "linewidth", 2 );
  xlabel ( "SNR-inj" );
  ylabel ( "log10(BSG)");
  %% ----------

  H_coverage = figure(); clf;
  perc = [PErecovery.f0_tau_percentile];
  perc = sort ( perc );
  Nn = length(perc);
  cdf = [1:Nn] / Nn;
  persistent covLU = [];
  if ( isempty ( covLU ) || any(size(covLU) != [2,Nn]) )
    DebugPrintf ( 1, "Computing new covLU[] ... ");
    covLU = zeros ( 2, Nn );
    for i = 1 : Nn
      [covLU(1, i), covLU(2, i)] = binomialConfidenceInterval ( Nn, i, 0.9 );
    endfor
    DebugPrintf ( 1, "done.\n");
  else
    DebugPrintf ( 1, "Re-using covLU[]\n");
  endif
  hold on;
  fill ( [ perc, fliplr(perc)], [covLU(1,:), fliplr(covLU(2,:))], "g", "facealpha", 0.5 );
  plot ( [0,1], [0,1], "--", "color", "black", "linewidth", 2 );
  stairs ( perc, cdf );
  legend ( "90% estimate", "exact" );

  grid on;
  xlabel ( "posterior percentile" );
  ylabel ( "measured coverage" );
  legend ( "location", "southeast" );

  return;

endfunction


function plotHists ( x, hists, varargin  = [] )

  numHists = length ( hists );
  assert ( length(x) == numHists );

  x_plot = p2sderr = p1sderr = meanerr = m1sderr = m2sderr = [];
  for l = 1 : numHists
    if ( isempty ( hists{l} ) ) continue; endif;
    h_l = hists{l}.hgrm;
    %% deal with singular (1-point or empty) histograms first
    count = histTotalCount ( h_l );
    if ( count == 0 ) continue; endif
    if ( all ( isinf ( histBins ( h_l, 1, "centre" ) ) ) ) continue; endif

    x_plot(end+1)  = x(l);

    if ( count == 1 )
      ii = find ( histProbs ( h_l ) );
      bin_ii = histBins( h_l, 1, "centre" )(ii);
      meanerr(end+1) = p2sderr(end+1) = p1sderr(end+1) = m1sderr(end+1) = m2sderr(end+1) = bin_ii;
    else
      ## calculate median and 2.5%, 25%, 75%, 97.5% quantiles
      p2sderr(end+1) = quantileFuncOfHist ( h_l, 0.025 );
      p1sderr(end+1) = quantileFuncOfHist ( h_l, 0.250 );
      meanerr(end+1) = quantileFuncOfHist ( h_l, 0.500 );
      m1sderr(end+1) = quantileFuncOfHist ( h_l, 0.750 );
      m2sderr(end+1) = quantileFuncOfHist ( h_l, 0.975 );
    endif
  endfor %% l = 1 : numHists

  ## do plots
  isHold = ishold();
  hold on;
  h1 = errorbar(x_plot, meanerr, meanerr - m1sderr, p1sderr - meanerr, "~");
  h2 = plot(x_plot, m2sderr, "--o", x_plot, p2sderr, "--o" );
  if ( !isHold ) hold; endif

  ## set plot properties
  h = [h1(:); h2(:)];
  set(h, "color", "black");
  set (h, "markersize", 2 );
  if length(varargin) > 0
    set(h, varargin{:});
  endif

  grid on;
endfunction
