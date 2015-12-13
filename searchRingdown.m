#!/usr/bin/octave -q

function searchRingdown ( varargin )
  global debugLevel = 1;

  uvar = parseOptions ( varargin,
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 400 },
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"tOffs", "real,strictpos,scalar", 0.43 },
                        {"prior_FreqRange", "real,strictpos,vector", [220,  280] },
                        {"prior_tauRange", "real,vector", [1e-3, 30e-3] },
                        {"prior_H", "real,strictpos,scalar", 1e-22},
                        {"showTSPlots", "bool", false }
                      );
  assert ( uvar.fMax > uvar.fMin );
  assert ( uvar.fMin < min(uvar.prior_FreqRange) );
  assert ( uvar.fMax > max(uvar.prior_FreqRange) );

  bname = sprintf ( "Ringdown-GPS%.0fs-f%.0fHz-%.0fHz-tau%.0fms-%.0fms-H%.2g",
                    uvar.tOffs, min(uvar.prior_FreqRange), max(uvar.prior_FreqRange),
                    min(uvar.prior_tauRange), max(uvar.prior_tauRange),
                    uvar.prior_H
                  );

  DebugPrintf ( 1, "Extracting timeseries ... ");
  [ts, tsW, tsOW, psd, IFO] = extractTS ( "fMin", uvar.fMin, "fMax", uvar.fMax, "tCenter", uvar.tCenter, "showPlots", uvar.showTSPlots, "simulate", "false", "lineSigma", 5, "lineWidth", 0.1 );
  DebugPrintf ( 1, "done.\n");
  Ndet = length(ts);

  %% ---------- estimate ~constant noise-floor level in band around ringdown frequencies ----------
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

  t0 = uvar.tCenter + uvar.tOffs;   	%% start-time of exponential ringdown-template
  Tmax = 5 * max(uvar.prior_tauRange);	%% max time range considered = 5 * tauMax

  %% ---------- templated search over {freq, tau} space ----------
  f0 = [min(uvar.prior_FreqRange):1:max(uvar.prior_FreqRange)];
  tau = [min(uvar.prior_tauRange):0.001:max(uvar.prior_tauRange)]';
  [ff, ttau] = meshgrid ( f0, tau );
  indsFreq = round ( (ff - fMin)/df );
  lap_s = 1./ ttau + I * 2*pi * ff;	%% laplace 'frequency'
  Gam_tot = ( ttau ./ Stot(indsFreq) + 1/uvar.prior_H^2 ).^(-1);
  BSG_tot = zeros ( size ( Gam_tot ) );
  match_tot = zeros ( size ( BSG_tot ) );
  DebugPrintf ( 1, "Searching ringdown %d templates ... ", length(ff(:)) );
  for X = 1:Ndet
    BSG{X} = match{X} = zeros ( size ( BSG_tot ) );

    inds   = find ( (tsOW{X}.ti >= t0) & (tsOW{X}.ti <= t0 + Tmax) );
    t_i    = tsOW{X}.ti(inds);
    Dt_i   = t_i - t0;
    xOW_i  = tsOW{X}.xi(inds);
    GamInv = ( ttau ./ psd{X}.Sn(indsFreq) + 1 / uvar.prior_H^2 );
    Gam{X} = 1 ./ GamInv;

    for l = 1 : length(lap_s(:))
      template_i = exp ( - Dt_i  * lap_s(l) );
      match{X}(l)   = 2 * dt * sum ( xOW_i .* template_i );
    endfor
    BSG{X} = compute_BSG ( Gam{X}, match{X}, uvar.prior_H );
    match_tot += match{X};

  endfor %% X = 1:Ndet
  BSG_tot = compute_BSG ( Gam_tot, match_tot, uvar.prior_H );
  DebugPrintf ( 1, "done.\n");

  DebugPrintf (2,  "<Gam> = %g, <Gam{1}> = %g, <Gam{2}> = %g\n", mean(Gam_tot(:)), mean(Gam{1}(:)), mean(Gam{2}(:)) );
  DebugPrintf (2,  "max|match| = %g, max|match{1}| = %g, max|match{2}| = %g\n", max(abs(match_tot(:))), max(abs(match{1}(:))), max(abs(match{2}(:))) );
  DebugPrintf (2,  "maxBSG = %g, maxBSG{1} = %g, maxBSG{2} = %g\n", max(BSG_tot(:)), max(BSG{1}(:)), max(BSG{2}(:)) );
  %% ---------- determine maximum-posterior estimates (MPE) ----------
  BSG_max = max(BSG_tot(:))
  l_MPE = (find ( BSG_tot(:) == BSG_max ))(1);
  f0_MPE = ff ( l_MPE )
  tau_MPE = ttau ( l_MPE )
  Gam_MPE = Gam_tot (l_MPE)
  F_MPE = match_tot(l_MPE)
  %% amplitude estimates:
  As_MPE = - Gam_MPE * imag ( F_MPE )
  Ac_MPE =   Gam_MPE * real ( F_MPE )
  A_MPE = sqrt ( As_MPE^2 + Ac_MPE^2 )
  phi0_MPE = atan2 ( -As_MPE, Ac_MPE )

  Dt = ts{1}.ti - t0;
  indsRingdown = find ( Dt >= 0 );
  tmpl_MPE = zeros ( size ( indsRingdown ) );
  Dt_pos = Dt(indsRingdown);
  tmpl_MPE = A_MPE * e.^(- Dt_pos / tau_MPE ) .* cos ( 2*pi * f0_MPE * Dt_pos + phi0_MPE );
  %% ---------- Plot results ----------
  %% Bayes factor + parameter-estimation posterior on {f,tau}
  figure(); clf;
  subplot ( 1+Ndet, 1, 1 );
  hold on;
  surf ( ff, ttau * 1e3, (BSG_tot) ); colorbar(); view(2); shading("interp");
  plot3 ( f0_MPE, tau_MPE * 1e3, 1.1*BSG_max, "marker", "x", "markersize", 3, "color", "white" );
  yrange = ylim(); xrange = xlim();
  xlabel ("Freq [Hz]"); ylabel ("$\\tau$ [ms]");
  tstring = sprintf ( "tOffs = +%.3f s: f in [%.0f, %.0f] Hz; BSG = %.2g", uvar.tOffs, uvar.fMin, uvar.fMax, mean(BSG_tot(:)) );
  title ( tstring );
  hold off;

  for X = 1:Ndet
    subplot ( 1+Ndet, 1, 1+X );
    hold on;
    surf ( ff, ttau * 1e3, (BSG{X}) ); colorbar(); view(2); shading("interp");
    plot3 ( f0_MPE, tau_MPE * 1e3, 1.1*BSG_max, "marker", "x", "markersize", 3, "color", "white" );
    yrange = ylim(); xrange = xlim();
    xlabel ("Freq [Hz]"); ylabel ("$\\tau$ [ms]");
    title ( tstring );
    hold off;
  endfor

  %% OW Time-series + MPE 'template'
  figure(); clf;
  %%subplot ( 2, 1, 2 );
  hold on;
  colors = { "red", "blue" };
  for X = 1:Ndet
    sleg = sprintf (";%s;", IFO{X} );
    plot ( ts{X}.ti - uvar.tCenter, ts{X}.xi, sleg, "linewidth", 2, "color", colors{X} );
  endfor
  plot ( ts{1}.ti(indsRingdown) - uvar.tCenter, tmpl_MPE, ";MPE;", "linewidth", 3, "color", "black" );
  yrange = ylim();
  line ( [ uvar.tOffs, uvar.tOffs], yrange, "linestyle", "--" );

  xlim ( [0.38, 0.46 ] );
  xlabel ("tOffs [s]");
  ylabel ("Strain/SX(f)");
  grid on;
  ylim ( yrange );
  hold off;

  return

endfunction

function BSG = compute_BSG ( Gam, match, H )
  Fstat = 0.5 * Gam .* abs ( match ).^2;
  BSG = Gam / H^2 .* e.^Fstat;
  return;
endfunction
