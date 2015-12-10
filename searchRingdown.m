#!/usr/bin/octave -q

function searchRingdown ( varargin )


  uvar = parseOptions ( varargin,
                        {"fMin", "real,strictpos,scalar", 200 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"tOffs", "real,strictpos,scalar", 0.43 },
                        {"prior_FreqRange", "real,strictpos,vector", [220,  280] },
                        {"prior_tauRange", "real,strictpos,vector", [2e-3, 100e-3] },
                        {"prior_H", "real,strictpos,scalar", 1e-22},
                        {"showPlots", "bool", false }
                      );
  assert ( uvar.fMax > uvar.fMin );

  [ts, tsW, tsOW, psd] = extractTS ( "fMin", uvar.fMin, "fMax", uvar.fMax, "tCenter", uvar.tCenter, "showPlots", uvar.showPlots, "simulate", "false" );
  Ndet = length(ts);

  SinvSum = 0;
  for X = 1:Ndet
    Sn{X} = mean ( psd{X}.Sn );
    SinvSum += 1/Sn{X};
    dt{X} = mean( diff ( ts{X}.ti ) );
  endfor
  Stot = 1 / ( 1/Ndet * SinvSum );

  assert ( max ( abs ( diff ( [dt{:}] ) )) < 1e-6 );
  dt = dt{1};

  t0 = uvar.tCenter + uvar.tOffs;   	%% start-time of exponential ringdown-template
  Tmax = 5 * max(uvar.prior_tauRange);	%% max time range considered = 5 * tauMax

  %% ---------- templated search over {freq, tau} space ----------
  f0 = [min(uvar.prior_FreqRange):1:max(uvar.prior_FreqRange)];
  tau = [min(uvar.prior_tauRange):0.001:max(uvar.prior_tauRange)]';
  [ff, ttau] = meshgrid ( f0, tau );
  lap_s = I * 2*pi * ff  + 1./ ttau;	%% laplace 'frequency'
  Gam_tot = ( ttau / Stot + 1/uvar.prior_H^2 ).^(-1);
  BSG_tot = zeros ( size ( Gam_tot ) );
  match_tot = zeros ( size ( BSG_tot ) );
  for X = 1:Ndet
    inds = find ( (tsOW{X}.ti >= t0) & (tsOW{X}.ti <= t0 + Tmax) );
    t_i   = tsOW{X}.ti(inds);
    xOW_i = tsOW{X}.xi(inds);
    Dt_i  = t_i - t0;
    BSG{X} = match{X} = zeros ( size ( BSG_tot ) );
    Gam{X} = ( ttau / Sn{X} + 1/uvar.prior_H^2 ).^(-1);
    for l = 1 : length(lap_s(:))
      template_i = exp ( - lap_s(l) * Dt_i );
      match{X}(l)   = dt * sum ( xOW_i .* template_i );
    endfor
    BSG{X} = compute_BSG ( Gam{X}, match{X}, uvar.prior_H );
    match_tot += match{X};

  endfor %% X = 1:Ndet
  BSG_tot = compute_BSG ( Gam_tot, match_tot, uvar.prior_H );

  figure(); clf; hold on;
  surf ( ff, ttau * 1e3, BSG_tot ); colorbar(); view(2); shading("interp");
  yrange = ylim(); xrange = xlim();
  xlabel ("Freq [Hz]"); ylabel ("$\\tau$ [ms]");
  tstring = sprintf ( "t0 = +%.3f s: f in [%.0f, %.0f] Hz; BSG = %.2g", uvar.tOffs, uvar.fMin, uvar.fMax, mean(BSG_tot(:)) );
  title ( tstring );
  grid on;
  fname = sprintf ( "PE-t%.3fs-freq%.0fHz-%.0fHz", uvar.tOffs, uvar.fMin, uvar.fMax );
  plot2pdf ( fname );
  return

endfunction

function BSG = compute_BSG ( Gam, match, H )
  Fstat = 0.5 * Gam .* abs ( match ).^2;
  BSG = Gam / H^2 .* e.^Fstat;
  return;
endfunction
