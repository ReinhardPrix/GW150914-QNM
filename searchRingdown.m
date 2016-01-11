#!/usr/bin/octave -q

function ret = searchRingdown ( varargin )
  global debugLevel = 1;

  uvar = parseOptions ( varargin,
                        {"ts", "cell" },	%% cell-array [over detectors]: normal, whitentend, and over-whitened timeseries
                        {"psd", "cell"},	%% cell-array [over detectors]: PSD estimate over frequency range, for each detector
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"tOffs", "real,scalar", 0.43 },
                        {"prior_f0Range", "real,strictpos,vector", [220,  270] },
                        {"prior_tauRange", "real,vector", [1e-3, 30e-3] },
                        {"prior_H", "real,strictpos,scalar", 4e-22},
                        {"plotResults", "bool", false }
                      );

  NoiseBand = 20;	%% use +-20Hz noise band around signal frequency for (non-OW) PSD estimate
  shiftL1 = 7.1e-3;	%% time-shift to apply to L1 data-stream

  assert ( min(uvar.psd{1}.fk) <= min(uvar.prior_f0Range) - NoiseBand );
  assert ( max(uvar.psd{1}.fk) >= max(uvar.prior_f0Range) + NoiseBand );

  bname = sprintf ( "Ringdown-GPS%.0fs-f%.0fHz-%.0fHz-tau%.1fms-%.1fms-H%.2g-tOffs%.4fs",
                    uvar.tCenter, min(uvar.prior_f0Range), max(uvar.prior_f0Range),
                    1e3 * min(uvar.prior_tauRange), 1e3 * max(uvar.prior_tauRange),
                    uvar.prior_H,
                    uvar.tOffs
                  );

  ts = uvar.ts; psd = uvar.psd;
  Ndet = length(ts);

  fk = psd{1}.fk;
  %% ---------- estimate ~constant noise-floor level in search-band of ringdown frequencies ----------
  indsSearch = find ( (fk > min(uvar.prior_f0Range)) & (fk < max(uvar.prior_f0Range)) );
  SinvSum = zeros ( size ( psd{1}.Sn ) );
  for X = 1:Ndet
    SinvSum += 1 ./ psd{X}.Sn;
    dt{X} = mean( diff ( ts{X}.ti ) );
    df{X} = mean ( diff ( psd{X}.fk ) );
    fMin{X} = min ( psd{X}.fk );
  endfor
  StotInv = (1/Ndet) * SinvSum;
  Stot = 1./ StotInv;

  assert ( max ( abs ( diff ( [dt{:}] ) )) < 1e-6 );
  assert ( max ( abs ( diff ( [df{:}] ) )) < 1e-6 );
  assert ( max ( abs ( diff ( [fMin{:}] ) )) < 1e-6 );
  dt = dt{1};
  df = df{1};
  fMin = fMin{1};

  Tmax = 3 * max(uvar.prior_tauRange);	%% max time range considered = 5 * tauMax

  %% ---------- templated search over {f0, tau} space ----------
  step_f0 = 0.1;
  step_tau = 0.0002;
  f0 = [min(uvar.prior_f0Range):step_f0:max(uvar.prior_f0Range)];
  tau = [min(uvar.prior_tauRange):step_tau:max(uvar.prior_tauRange)]';
  Nf0 = length(f0);
  Ntau = length(tau);
  [ff0, ttau] = meshgrid ( f0, tau );
  lap_s = 1./ ttau + I * 2*pi * ff0;	%% laplace 'frequency'
  Ntempl = length ( lap_s(:) );
  log10BSG_tot = match_tot = Gam_tot = zeros ( size ( lap_s ) );
  Is = Ic = Isc = det_gam = normBSG = zeros ( size ( lap_s ) );
  gam_ss = gam_cc = gam_sc = zeros ( size ( lap_s ) );
  DebugPrintf ( 1, "Searching ringdown %d templates ... ", length(ff0(:)) );

  %% ----- prepare time-series stretch to analyze
  t0 = uvar.tCenter + uvar.tOffs;   	%% start-time of exponential ringdown-template

  for X = 1:Ndet
    log10BSG{X} = match{X} = Gam{X} = zeros ( size ( lap_s ) );

    %% prepare per-detector timeseries for matching: adapt L1 data to be phase-coherent with H1
    if ( strcmp ( ts{X}.IFO, "L1" ) )
      shiftL1_eff = round ( shiftL1 / dt ) * dt;
      ts{X}.ti += shiftL1_eff;	%% make sure we shift by integer number of time-samples
      ts{X}.xi   *= -1;
      ts{X}.xiW  *= -1;
      ts{X}.xiOW *= -1;
    endif

    inds_match = find ( (ts{X}.ti >= t0) & (ts{X}.ti <= t0 + Tmax) );
    Dt_i{X} = ts{X}.ti ( inds_match ) - t0;	%% we made sure time-steps should agree (to within 2e-7s)
    xOW_i{X} = ts{X}.xiOW ( inds_match );
  endfor %% X

  %% ----- determine ~const noise-estimate in f0 +- 20Hz Band around each signal frequency f0
  inds_f0   = round ( (f0 - fMin)/df );
  inds_Band = round( NoiseBand / df);
  offs = [-inds_Band : inds_Band ];
  [xx, yy] = meshgrid ( inds_f0, offs );
  inds_mat = xx + yy;
  Stot_k  = mean ( Stot ( inds_mat ), 1 );
  Stot_mat = meshgrid ( Stot_k, ones(1,Ntau) );
  Gam_tot = ( ttau ./ Stot_mat  + 1 / uvar.prior_H^2 ).^(-1);
  for X = 1:Ndet
    Sn_k{X} = mean ( psd{X}.Sn ( inds_mat ), 1 );
    Sn_mat{X} = meshgrid ( Sn_k{X}, ones(1,Ntau) );
    Gam{X} = ( ttau ./ Sn_mat{X} + 1 / uvar.prior_H^2 ).^(-1);
  endfor

  %% ---------- search parameter-space in {f0, tau} and compute matched-filter in each template ----------
  for l = 1 : Ntempl	%% loop over all templates
    %% ----- (complex) time-domain template
    hExp_i = exp ( - Dt_i{1}  * lap_s(l) );
    for X = 1:Ndet
      match{X}(l) = 2 * dt * sum ( xOW_i{X} .* hExp_i );
    endfor %% X
  endfor
  for X = 1:Ndet
    match_tot   += match{X};
    x_s = - imag ( match_tot );
    x_c =   real ( match_tot );
  endfor

  for l = 1 : Ntempl
    %% ----- whitened frequency-domain template basis functions ----------
    denom_k = ( 1 + I * 4*pi * fk * ttau(l) - 4*pi^2 * ( fk.^2 - ff0(l)^2 ) * ttau(l)^2 ) .* sqrt ( Stot );
    hsFT_k  = ttau(l) * ( 2*pi * ff0(l) * ttau(l) ) ./ denom_k;
    hcFT_k  = ttau(l) * ( 1 + I * 2*pi*fk * ttau(l) ) ./ denom_k;
    %% ----- compute M-matrix from template self-match integrals in frequency-domain ----------
    Is(l)  = 4 * Ndet * df * sum ( abs(hsFT_k).^2 );
    Ic(l)  = 4 * Ndet * df * sum ( abs(hcFT_k).^2 );
    Isc(l) = 4 * Ndet * df * real ( sum ( hsFT_k .* conj(hcFT_k) ) );
  endfor %% l
  detM = Is .* Ic - Isc.^2;
  trM  = Is + Ic;
  Hm2 = uvar.prior_H^(-2);
  det_gamInv = ( Is + Hm2 ) .* ( Ic + Hm2 ) - Isc.^2;
  det_gam = 1./ det_gamInv;
  normBSG = sqrt(det_gam) * Hm2;
  gam_ss  = det_gam .* ( Ic + Hm2 );
  gam_cc  = det_gam .* ( Is + Hm2 );
  gam_sc  = det_gam .* ( -Isc );

  x_gam_x = gam_ss .* x_s.^2 + 2 * gam_sc .* x_s .* x_c + gam_cc .* x_c.^2;
  BSGalt = normBSG .* exp ( 0.5 .* x_gam_x );
  for X = 1:Ndet
    log10BSG{X} = compute_log10BSG ( Gam{X}, match{X}, uvar.prior_H );
    BSG{X} = 10.^(log10BSG{X});
  endfor
  log10BSG_tot = compute_log10BSG ( Gam_tot, match_tot, uvar.prior_H );
  BSG_tot = 10.^log10BSG_tot;


  posterior = struct ( "f0", ff0, "tau", ttau, "BSG", BSG_tot );
  DebugPrintf ( 1, "done.\n");

  DebugPrintf (2,  "<Gam> = %g, <Gam{1}> = %g, <Gam{2}> = %g\n", mean(Gam_tot(:)), mean(Gam{1}(:)), mean(Gam{2}(:)) );
  DebugPrintf (2,  "max|match| = %g, max|match{1}| = %g, max|match{2}| = %g\n", max(abs(match_tot(:))), max(abs(match{1}(:))), max(abs(match{2}(:))) );
  DebugPrintf (2,  "max[log10BSG] = %g, max[log10BSG{1}] = %g, max[log10BSG{2}] = %g\n", max(log10BSG_tot(:)), max(log10BSG{1}(:)), max(log10BSG{2}(:)) );

  %% ---------- compute marginalized posteriors on {f,tau} ----------
  %% renormalize to avoid overflow
  log10BSG_renorm = log10BSG_tot - max(log10BSG_tot(:));	%% == divide BSGL by max(BSGL) ==> max now at 1
  BSG_renorm = 10.^log10BSG_renorm;
  posterior_f0   = sum ( BSG_renorm, 1 );
  norm_f0 = step_f0 * sum ( posterior_f0 );
  posterior_f0 /= norm_f0;

  posterior_tau = sum ( BSG_renorm, 2 );
  norm_tau = step_tau * sum ( posterior_tau);
  posterior_tau /= norm_tau;

  %% find 90% credible intervals
  confidence = 0.90;
  try
    f0_est  = credibleInterval ( f0, posterior_f0, confidence );
    tau_est = credibleInterval ( tau, posterior_tau, confidence );
    failed_intervals = false;
  catch
    failed_intervals = true;
  end_try_catch

  %% ---------- estimate SNR for all templates ----------
  AVec_est = Gam_tot .* match_tot;
  A_est = abs ( AVec_est );
  phi0_est = arg ( AVec_est );
  SNR_est = A_est .* sqrt ( Ndet * ttau ./ (2 * Stot_mat ) );

  %% ---------- determine maximum-posterior estimates (MPE) ----------
  BSG_max = max ( BSG_tot(:) );
  l_MPE = ( find ( BSG_tot(:) == BSG_max ) )(1);
  Stot_MPE = Stot_mat ( l_MPE );
  f0_MPE = ff0 ( l_MPE );
  tau_MPE = ttau ( l_MPE );
  Gam_MPE = Gam_tot ( l_MPE );
  F_MPE = match_tot ( l_MPE );

  A_MPE = A_est(l_MPE);
  phi0_MPE = phi0_est(l_MPE);
  SNR_MPE = SNR_est(l_MPE);

  %% ---------- Plot results ----------

  %% ----- Fig 1: Bayes factor / posterior over {f0,tau} ----------
  if ( uvar.plotResults )
    figure(1); clf;
    subplot ( 2, 2, 1 );
    hold on;
    colormap ("jet");
    surf ( ff0, ttau * 1e3, BSG_tot ); colorbar("location", "NorthOutside"); view(2); shading("interp");
    plot3 ( f0_MPE, tau_MPE * 1e3, 1.1*BSG_max, "marker", "o", "markersize", 3, "color", "white" );
    yrange = ylim();
    xlabel ("f0 [Hz]"); ylabel ("tau [ms]");
    hold off;

    subplot ( 2, 2, 3 );
    plot ( f0, posterior_f0, "linewidth", 2 );
    grid on;
    yrange = ylim();
    if ( !failed_intervals )
      line ( [f0_est.MPE, f0_est.MPE], yrange );
      line ( [f0_est.lower, f0_est.upper], [f0_est.pIso, f0_est.pIso] );
    endif
    ylim ( [0, max(yrange)] );
    xlabel ("f0 [Hz]");
    ylabel ("pdf(f0)");

    subplot ( 2, 2, 2 );
    plot ( tau * 1e3, posterior_tau, "linewidth", 2 );
    yrange = ylim();
    if ( !failed_intervals )
      line ( [tau_est.MPE, tau_est.MPE]*1e3, yrange );
      line ( [tau_est.lower, tau_est.upper]* 1e3, [tau_est.pIso, tau_est.pIso] );
    endif
    ylim ( [0, max(yrange)] );
    grid on;
    xlabel ("tau [ms]");
    ylabel ("pdf(tau)");

    subplot ( 2, 2, 4 );
    hold on;
    colors = { "red", "blue" };
    for X = 1:Ndet
      sleg = sprintf (";%s;", ts{X}.IFO );
      plot ( ts{X}.ti - uvar.tCenter, ts{X}.xiOW * Stot_MPE, sleg, "linewidth", 2, "color", colors{X} );
    endfor
    Dt = ts{1}.ti - t0;
    indsRingdown = find ( Dt >= 0 );
    Dt_pos = Dt ( indsRingdown );
    tmpl_MPE = A_MPE * e.^(- Dt_pos / tau_MPE ) .* cos ( 2*pi * f0_MPE * Dt_pos + phi0_MPE );
    plot ( ts{1}.ti(indsRingdown) - uvar.tCenter, tmpl_MPE, ";MPE;", "linewidth", 3, "color", "black" );
    legend ( "location", "NorthEast");
    yrange = ylim();
    line ( [ uvar.tOffs, uvar.tOffs], yrange, "linestyle", "-", "linewidth", 3 );
    xlim ( [uvar.tOffs - 0.3 * Tmax, uvar.tOffs + Tmax ] );
    ylim ( yrange );
    xlabel ( sprintf ( "%.0f + tOffs [s]", uvar.tCenter) );
    tOffs_text = sprintf ( "tOffs = %.4f s", uvar.tOffs );
    x0 = uvar.tOffs + 0.02 * abs(diff(xlim()));
    y0 = max(yrange) - 0.2*abs(diff(yrange));
    text ( x0, y0, tOffs_text );
    text ( min(xlim()) - 0.2 * abs(diff(xlim())), 0, "h(t)" );
    grid on;
    hold off;

    figure(); clf;
    hold on;
    colormap ("jet");
    surf ( ff0, ttau * 1e3, BSGalt ); colorbar("location", "NorthOutside"); view(2); shading("interp");
    yrange = ylim();
    xlabel ("f0 [Hz]"); ylabel ("tau [ms]");
    hold off;

  endif %% plotResults

  %% summarize numerical outcomes on stdout
  BSG_mean = mean ( BSG_tot(:) );
  BSGalt_mean = mean ( BSGalt(:) );

  DebugPrintf (1, "tGPS    = %.0f + %f s\n", uvar.tCenter, uvar.tOffs );
  DebugPrintf (1, "<BSG>   = %.2g\n", BSG_mean );
  DebugPrintf (1, "<BSGalt>= %.2g\n", BSGalt_mean );
  DebugPrintf (1, "f0_est  = { %.1f, %.1f, %1.f } Hz\n", f0_est.lower, f0_est.MPE, f0_est.upper );
  DebugPrintf (1, "tau_est = { %.1f, %.1f, %1.f } ms\n", 1e3 * tau_est.lower, 1e3 * tau_est.MPE, 1e3 * tau_est.upper );
  DebugPrintf (1, "A_MPE   = %.2g\n", A_MPE );
  DebugPrintf (1, "phi0_MPE= %.2g\n", phi0_MPE );
  DebugPrintf (1, "sqrt(S)(f_MPE)= %.2g\n", sqrt(Stot_MPE) );
  DebugPrintf (1, "SNR_MPE = %.2g\n", SNR_MPE );
  DebugPrintf (1, "sanity check: BSG_s = Lr_s ~ e^(SNR^2/2) = %.2g\n", exp ( SNR_MPE^2 / 2 ) );


  ret = struct ( "bname", bname, ...
                 "tGPS", uvar.tCenter + uvar.tOffs, ...
                 "BSG_mean", BSG_mean, ...
                 "A_MPE", A_MPE, ...
                 "phi0_MPE", phi0_MPE, ...
                 "f0_est", f0_est, ...
                 "tau_est", tau_est, ...
                 "SNR_MPE", SNR_MPE, ...
                 "sqrtS_MPE", sqrt(Stot_MPE), ...
                 "posterior", posterior, ...
                 "SNR", SNR_est
               );

  return

endfunction

function log10BSG = compute_log10BSG ( Gam, match, H )
  Fstat = 0.5 * Gam .* abs ( match ).^2;
  log10BSG = log10(Gam) - 2 * log10(H) + Fstat * log10(e);
  return;
endfunction


%% return MP estimate 'x_est' with fields 'MPE', 'lower', 'upper', and 'pIso'
%% ie: maximum-posterior estimate x_est.MPE
%% 'confidence'-credible interval [x_est.lower, x_est.upper], and
%% corresponding iso-posterior value x_est.pIso
function x_est = credibleInterval ( x, posterior_x, confidence = 0.9 )

  assert ( size ( x ) == size ( posterior_x ) );
  dx = mean ( diff ( x ) );
  [x_pIso, delta, INFO, OUTPUT] = fzero ( @(x_pIso)  dx * sum ( posterior_x ( find ( posterior_x >= x_pIso ) ) ) - confidence, ...
                                          [ min(posterior_x), max(posterior_x) ], ...
                                          optimset ( "TolX", 1e-4 )
                                        );
  try
    assert ( INFO == 1 );
  catch
    delta
    INFO
    OUTPUT
    error ("fzero() failed\n");
  end_try_catch
  x_MP = x ( find ( posterior_x == max(posterior_x) ) );

  inds0 = find ( posterior_x >= x_pIso );
  i_min = min(inds0);
  i_max = max(inds0);
  %% check if these are respective closest to p_iso
  if ( (i_min > 1) && abs(posterior_x(i_min-1)-x_pIso) < abs(posterior_x(i_min)-x_pIso) )
    i_min --;
  endif
  if ( (i_max < length(posterior_x)) && abs(posterior_x(i_max+1)-x_pIso) < abs(posterior_x(i_max)-x_pIso) )
    i_max ++;
  endif

  x_est = struct ( "MPE", x_MP, "lower", x(i_min), "upper", x(i_max), "pIso", x_pIso );
  return;
endfunction
