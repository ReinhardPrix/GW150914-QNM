#!/usr/bin/octave -q

function ret = searchRingdown ( varargin )
  global debugLevel = 1;

  uvar = parseOptions ( varargin,
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 400 },
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"tOffs", "real,strictpos,scalar", 0.43 },
                        {"prior_FreqRange", "real,strictpos,vector", [220,  280] },
                        {"prior_tauRange", "real,vector", [1e-3, 30e-3] },
                        {"prior_H", "real,strictpos,scalar", 1e-22},
                        {"extraLabel", "char,vector", ""},
                        {"showTSPlots", "bool", false },
                        {"plotResults", "bool", false },
                        {"useOW", "bool", true}	%% use overwhitening to handle matched-filtering noise estimate, alt: estimate ~const noise floor
                      );
  assert ( uvar.fMax > uvar.fMin );
  assert ( uvar.fMin < min(uvar.prior_FreqRange) );
  assert ( uvar.fMax > max(uvar.prior_FreqRange) );

  bname = sprintf ( "Ringdown-GPS%.0fs-f%.0fHz-%.0fHz-tau%.0fms-%.0fms-H%.2g-%s-%s",
                    uvar.tCenter, min(uvar.prior_FreqRange), max(uvar.prior_FreqRange),
                    min(uvar.prior_tauRange), max(uvar.prior_tauRange),
                    uvar.prior_H,
                    ifelse ( uvar.useOW, "OW", "constS" ),
                    uvar.extraLabel
                  );

  DebugPrintf ( 1, "Extracting timeseries ... ");
  [ts, tsW, tsOW, psd, IFO] = extractTS ( "fMin", uvar.fMin, "fMax", uvar.fMax, "tCenter", uvar.tCenter, "showPlots", uvar.showTSPlots, "lineSigma", 5, "lineWidth", 0.1, "RngMedWindow", 300 );
  DebugPrintf ( 1, "done.\n");
  Ndet = length(ts);

  %% ---------- estimate ~constant noise-floor level in search-band of ringdown frequencies ----------
  indsSearch = find ( (psd{1}.fk > min(uvar.prior_FreqRange)) & (psd{1}.fk < max(uvar.prior_FreqRange)) );
  SinvSum = zeros ( size ( psd{1}.Sn ) );
  for X = 1:Ndet
    SinvSum += 1 ./ psd{X}.Sn;
    dt{X} = mean( diff ( tsOW{X}.ti ) );
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

  %% ---------- templated search over {freq, tau} space ----------
  step_f0 = 0.1;
  step_tau = 0.0002;
  f0 = [min(uvar.prior_FreqRange):step_f0:max(uvar.prior_FreqRange)];
  tau = [min(uvar.prior_tauRange):step_tau:max(uvar.prior_tauRange)]';
  Nf0 = length(f0);
  Ntau = length(tau);
  [ff, ttau] = meshgrid ( f0, tau );
  lap_s = 1./ ttau + I * 2*pi * ff;	%% laplace 'frequency'
  log10BSG_tot = match_tot = Gam_tot = zeros ( size ( lap_s ) );
  DebugPrintf ( 1, "Searching ringdown %d templates ... ", length(ff(:)) );

  %% ----- prepare time-series stretch to analyze
  t0 = uvar.tCenter + uvar.tOffs;   	%% start-time of exponential ringdown-template
  inds = find ( (tsOW{1}.ti >= t0) & (tsOW{1}.ti <= t0 + Tmax) );
  t_i  = ts{1}.ti(inds);
  Dt_i = t_i - t0;

  for X = 1:Ndet
    log10BSG{X} = match{X} = Gam{X} = zeros ( size ( lap_s ) );
    xOW_i{X} = tsOW{X}.xi(inds);
    x_i{X}   = ts{X}.xi(inds);
  endfor %% X

  %% ----- determine ~const noise-estimate in f0 +- 20Hz Band around each signal frequency f0
  inds_f0   = round ( (f0 - fMin)/df );
  inds_Band = round(20 / df);
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
  Ntempl = length ( lap_s(:) );
  for l = 1 : Ntempl	%% loop over all templates
    template_l = exp ( - Dt_i  * lap_s(l) );
    for X = 1:Ndet
      if ( uvar.useOW )	%% use overwhitened time-series for match
        match{X}(l) = 2 * dt * sum ( xOW_i{X} .* template_l );
      else	%% use ~const noise-floor estimate in signal +-20Hz band
        match{X}(l) = 2 * dt * sum ( x_i{X} .* template_l ) / Sn_mat{X}(l);
      endif

    endfor %% X
  endfor %% l

  for X = 1:Ndet
    match_tot   += match{X};
    log10BSG{X} = compute_log10BSG ( Gam{X}, match{X}, uvar.prior_H );
    BSG{X} = 10.^(log10BSG{X});
  endfor
  log10BSG_tot = compute_log10BSG ( Gam_tot, match_tot, uvar.prior_H );
  BSG_tot = 10.^log10BSG_tot;

  DebugPrintf ( 1, "done.\n");

  DebugPrintf (2,  "<Gam> = %g, <Gam{1}> = %g, <Gam{2}> = %g\n", mean(Gam_tot(:)), mean(Gam{1}(:)), mean(Gam{2}(:)) );
  DebugPrintf (2,  "max|match| = %g, max|match{1}| = %g, max|match{2}| = %g\n", max(abs(match_tot(:))), max(abs(match{1}(:))), max(abs(match{2}(:))) );
  DebugPrintf (2,  "max[log10BSG] = %g, max[log10BSG{1}] = %g, max[log10BSG{2}] = %g\n", max(log10BSG_tot(:)), max(log10BSG{1}(:)), max(log10BSG{2}(:)) );
  %% ---------- determine maximum-posterior estimates (MPE) ----------
  BSG_max = max(BSG_tot(:));
  l_MPE = (find ( BSG_tot(:) == BSG_max ))(1);
  Stot_MPE = Stot_mat( l_MPE );
  f0_MPE = ff ( l_MPE );
  tau_MPE = ttau ( l_MPE );
  Gam_MPE = Gam_tot (l_MPE);
  F_MPE = match_tot(l_MPE);
  %% amplitude estimates:
  As_MPE = - Gam_MPE * imag ( F_MPE );
  Ac_MPE =   Gam_MPE * real ( F_MPE );
  A_MPE = sqrt ( As_MPE^2 + Ac_MPE^2 );
  phi0_MPE = atan2 ( -As_MPE, Ac_MPE );

  Dt = ts{1}.ti - t0;
  indsRingdown = find ( Dt >= 0 );
  tmpl_MPE = zeros ( size ( indsRingdown ) );
  Dt_pos = Dt(indsRingdown);
  tmpl_MPE = A_MPE * e.^(- Dt_pos / tau_MPE ) .* cos ( 2*pi * f0_MPE * Dt_pos + phi0_MPE );

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
  f0_est  = credibleInterval ( f0, posterior_f0, confidence );
  tau_est = credibleInterval ( tau, posterior_tau, confidence );
  %% ---------- Plot results ----------

  %% ----- Fig 1: Bayes factor / posterior over {f0,tau} ----------
  if ( uvar.plotResults )
    figure(); clf;
    subplot ( 2, 2, 1 );
    hold on;
    colormap ("jet");
    surf ( ff, ttau * 1e3, BSG_tot ); colorbar("location", "NorthOutside"); view(2); shading("interp");
    plot3 ( f0_MPE, tau_MPE * 1e3, 1.1*BSG_max, "marker", "o", "markersize", 3, "color", "white" );
    yrange = ylim();
    xlabel ("Freq [Hz]"); ylabel ("tau [ms]");
    hold off;

    subplot ( 2, 2, 3 );
    plot ( f0, posterior_f0, "linewidth", 2 );
    grid on;
    yrange = ylim();
    line ( [f0_est.MPE, f0_est.MPE], yrange );
    line ( [f0_est.lower, f0_est.upper], [f0_est.pIso, f0_est.pIso] );
    ylim ( yrange );
    xlabel ("Freq [Hz]");
    ylabel ("pdf(Freq)");

    subplot ( 2, 2, 2 );
    plot ( tau * 1e3, posterior_tau, "linewidth", 2 );
    yrange = ylim();
    line ( [tau_est.MPE, tau_est.MPE]*1e3, yrange );
    line ( [tau_est.lower, tau_est.upper]* 1e3, [tau_est.pIso, tau_est.pIso] );
    ylim ( yrange );
    grid on;
    xlabel ("tau [ms]");
    ylabel ("pdf(tau)");

    subplot ( 2, 2, 4 );
    hold on;
    colors = { "red", "blue" };
    for X = 1:Ndet
      sleg = sprintf (";%s;", IFO{X} );
      plot ( ts{X}.ti - uvar.tCenter, ts{X}.xi, sleg, "linewidth", 2, "color", colors{X} );
    endfor
    plot ( ts{1}.ti(indsRingdown) - uvar.tCenter, tmpl_MPE, ";MPE;", "linewidth", 3, "color", "black" );
    legend ( "location", "NorthWest");
    yrange = ylim();
    line ( [ uvar.tOffs, uvar.tOffs], yrange, "linestyle", "-", "linewidth", 2 );
    xlim ( [0.38, 0.46 ] );
    ylim ( yrange );
    xlabel ("tOffs [s]");
    tOffs_text = sprintf ( "tOffs = %.3f s", uvar.tOffs );
    x0 = uvar.tOffs + 0.02 * abs(diff(xlim()));
    y0 = max(yrange) - 0.2*abs(diff(yrange));
    text ( x0, y0, tOffs_text );
    text ( min(xlim()) - 0.2 * abs(diff(xlim())), 0, "h(t)" );
    grid on;
    hold off;

    %% ---------- output final results summary in input-specific collection file ----------
    fname = sprintf ( "%s.png", bname );
    ezprint ( fname, "width", 512 );

  endif %% plotResults

  %% summarize numerical outcomes on stdout
  BSG_mean = mean ( BSG_tot(:) );

  DebugPrintf (1, "tGPS    = %.0f + %f s\n", uvar.tCenter, uvar.tOffs );
  DebugPrintf (1, "<BSG>   = %.2g\n", BSG_mean );
  DebugPrintf (1, "f0_est  = { %.1f, %.1f, %1.f } Hz\n", f0_est.lower, f0_est.MPE, f0_est.upper );
  DebugPrintf (1, "tau_est = { %.1f, %.1f, %1.f } ms\n", 1e3 * tau_est.lower, 1e3 * tau_est.MPE, 1e3 * tau_est.upper );
  DebugPrintf (1, "A_MPE   = %.2g\n", A_MPE );
  DebugPrintf (1, "phi0_MPE= %.2g\n", phi0_MPE );
  DebugPrintf (1, "sqrt(S)(f_MPE)= %.2g\n", sqrt(Stot_MPE) );
  SNR_MPE = A_MPE * sqrt ( Ndet * tau_est.MPE / (2 * Stot_MPE ) );
  DebugPrintf (1, "SNR_MPE = %.2g\n", SNR_MPE );
  DebugPrintf (1, "sanity check: BSG_s = Lr_s ~ e^(SNR^2/2) = %.2g\n", exp ( SNR_MPE^2 / 2 ) );


  ret = struct ( "tGPS", uvar.tCenter + uvar.tOffs, ...
                 "BSG_mean", BSG_mean, ...
                 "A_MPE", A_MPE, ...
                 "phi0_MPE", phi0_MPE, ...
                 "f0_est", f0_est, ...
                 "tau_est", tau_est, ...
                 "SNR_MPE", SNR_MPE, ...
                 "sqrtS_MPE", sqrt(Stot_MPE)
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
