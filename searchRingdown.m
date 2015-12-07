#!/usr/bin/octave -q

function searchRingdown ( varargin )


  uvar = parseOptions ( varargin,
                        {"fMin", "real,strictpos,scalar", 200 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"tOffs", "real,strictpos,scalar", 0.43 },
                        {"extraLabel", "char,vector", ""},
                        {"prior_FreqRange", "real,strictpos,vector", [220,  280] },
                        {"prior_tauRange", "real,strictpos,vector", [2e-3, 100e-3] }
                      );
  assert ( uvar.fMax > uvar.fMin );

  IFOs = { "H1"; "L1" };
  for X = 1:length(IFOs)
    fnames{X} = sprintf ( "TS-%s-freq%.0fHz-%.0fHz-20s%s.dat", IFOs{X}, uvar.fMin, uvar.fMax, uvar.extraLabel );
    dat = load( fnames{X} );
    tS{X}.ti = dat(:,1);
    tS{X}.xi = dat(:,2);
  endfor

  Sn_H1 = estimatePSD ( tS{1}.ti, tS{1}.xi, min(uvar.prior_FreqRange), max(uvar.prior_FreqRange) );
  Sn_L1 = estimatePSD ( tS{2}.ti, tS{2}.xi, min(uvar.prior_FreqRange), max(uvar.prior_FreqRange) );
  Sn = harmmean ( [Sn_H1, Sn_L1] );
  printf ( "sqrtSn[H1] = %.2g, sqrtSn[L1] = %.2g, sqrtSn = %.2g\n", sqrt(Sn_H1), sqrt(Sn_L1), sqrt(Sn) );

  dt_H1 = mean( diff ( tS{1}.ti ) );
  dt_L1 = mean( diff ( tS{2}.ti ) );
  assert ( dt_L1, dt_H1, 1e-6 );
  dt = dt_H1;

  tEvent = 1126259462;

  H = 1e-22;	%% prior scale for signal strength
  T = 0.25;	%% max time range considered

  %% extract time-window of interest, [t0, t0+T], bring both IFOs into phase-coherence (L1+7.3ms and x(-1)), and sum
  tShift0_L1 = 7.3e-3;	%% L1 arrival time 7.3ms earlier: shift forward
  iShift_L1 = round ( tShift0_L1 / dt );
  tShift_L1 = iShift_L1 * dt; %% actual time-shift
  assert ( tShift0_L1, tShift_L1, 1e-4 );	%% make sure we got close enough, tolerance 0.1ms

  %% start-time of exponential decay-window
  t0 = tEvent + uvar.tOffs;

  inds = find ( (tS{1}.ti >= t0) & (tS{1}.ti <= t0 + T) );
  ti = tS{1}.ti(inds);	%% time 'reference' is arrival time at H1

  xi_H1 = tS{1}.xi(inds);
  xi_L1 = (-1)*tS{2}.xi(inds - iShift_L1); %% shift this into the future, means use earlier indices as later ones

  %% combine coherent data-streams, weighted by respective noise-floors
  yi = xi_H1 / Sn_H1 + xi_L1 / Sn_L1;

  f0 = [min(uvar.prior_FreqRange):1:max(uvar.prior_FreqRange)];
  tau = [min(uvar.prior_tauRange):0.001:max(uvar.prior_tauRange)]';
  [ff, ttau] = meshgrid ( f0, tau );
  %% Laplace transform for frequency 'f0' and damping factor 'Q': s = 2pi f0 * (i + 1/Q)
  Q = pi * ff .* ttau;

  s = 2*pi * ff .* ( I + 1./(2*Q) );
  BF = zeros ( size ( s ) );

  for l = 1 : length(s(:))
    template = exp ( -s(l) * (ti - t0) );
    match = 2 * dt * sum ( yi .* template );
    Gam = ( ttau(l) / (2*Sn) + 1/H^2 )^(-1);
    Fstat = 0.5 * Gam * abs ( match )^2;
    BF(l) = (1 + H^2 * ttau(l) / ( 2 * Sn ) )^(-1) * e^Fstat;
  endfor

  BFmean = mean ( BF(:) );
  BFmax  = max(max(BF));
  [imax, jmax] = find ( BF == BFmax );
  f0MPE = ff ( imax, jmax );
  tauMPE = ttau ( imax, jmax );

  figure(); clf; hold on;
  surf ( ff, ttau * 1e3, BF ); colorbar(); view(2); shading("interp");
  yrange = ylim(); xrange = xlim();
  line ( [f0MPE, f0MPE], yrange, "linestyle", "--" );
  line ( xrange, 1e3*[tauMPE, tauMPE], "linestyle", "--" );
  xlabel ("Freq [Hz]"); ylabel ("$\\tau$ [ms]");
  tstring = sprintf ( "t0 = +%.3f s: f in [%.0f, %.0f] Hz; BSG = %.2g\n", uvar.tOffs, uvar.fMin, uvar.fMax, BFmean );
  title ( tstring );
  grid on;
  fname = sprintf ( "PE-t%.3fs-freq%.0fHz-%.0fHz%s", uvar.tOffs, uvar.fMin, uvar.fMax, uvar.extraLabel );
  plot2pdf ( fname );
  return
endfunction

%% estimate single-sided noise PSD 'Sn':
function Sn = estimatePSD ( ti, xi, fMin, fMax )
  ft = FourierTransform ( ti, xi' );
  dt = mean ( diff ( ti ) );
  T = max(ti) - min(ti) + dt;
  inds = find ( (ft.fk >= fMin) & (ft.fk <= fMax) );
  periodo = abs(ft.xk(inds)).^2;
  %% Wiener-Khintchine theorm
  %% Sn_mean = 1/T * mean ( periodo );
  try
    lal;
    Sn  = 2/T * median ( periodo ) / XLALRngMedBias ( length(periodo));
  catch
    error ("Need SWIG Lalsuite installation for XLALRngMedBias()\n");
  end_try_catch
  return;
endfunction
