%% estimate single-sided noise PSD 'psd(f)', and return whitened and over-whitened TS
%% including an automated line-nuking algorithm: 'lines' are identified as > lineSigma deviations in power
function [ psd, tsOut ] = whitenTS ( varargin )
  global debugLevel = 1;

  uvar = parseOptions ( varargin,
                        {"tsIn", "struct" },
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"lineSigma", "real,positive,scalar", 5},	%% sigma deviations to indentify 'lines' in spectrum
                        {"lineWidth", "real,positive,scalar", 0.1},	%% +- width in Hz to zero around 'lines'
                        {"RngMedWindow", "real,positive,scalar", 300 },  %% window size to use for rngmed-based PSD estimation
                        {"plotSpectrum", "bool", false }
                      );

  ti = uvar.tsIn.ti;
  tStart = min ( ti );
  dt = mean ( diff ( ti ) );
  fSamp = 1/dt;
  T = max(ti) - tStart + dt;
  ft = FourierTransform ( uvar.tsIn.ti, uvar.tsIn.xi );

  sideband = uvar.RngMedWindow / (2*T); 		%% extra frequency side-band for median PSD estimates
  indsWide = find ( (ft.fk >=  uvar.fMin - sideband) & (ft.fk <=  uvar.fMax + sideband) );
  indsUse  = find ( (ft.fk(indsWide) >= uvar.fMin) & (ft.fk(indsWide) <=  uvar.fMax ) );

  %% ---------- estimate PSD using running-median over frequency ----------
  %% Wiener-Khintchine theorm
  %% Sn_double = 1/T * mean ( periodo );
  lal;
  rngmedbias = XLALRngMedBias ( uvar.RngMedWindow );

  periodo = abs ( ft.xk(indsWide) ).^2;
  Sn_wide = 2/T * rngmed ( periodo, uvar.RngMedWindow ) / rngmedbias;	%% single-sided PSD

  psd.fk = ft.fk ( indsWide ( indsUse ) );
  psd.Sn = Sn_wide ( indsUse );
  psd.IFO   = uvar.tsIn.IFO;

  %% ---------- automatically identify and nuke 'lines' ----------
  fk_wide   = ft.fk ( indsWide );
  xk_wide   = ft.xk ( indsWide );

  yk = xk_wide  ./ sqrt( T * Sn_wide / 2 );   %% normalized SFT:
  Pk_re = abs ( real ( yk ) );
  Pk_im = abs ( imag ( yk ) );
  inds_lines = find ( (Pk_re > uvar.lineSigma) | (Pk_im > uvar.lineSigma) );
  iNukeWidth = round ( uvar.lineWidth * T );
  sides = -iNukeWidth : iNukeWidth;
  [xx, yy] = meshgrid ( inds_lines, sides );
  nuke_lines = xx + yy;
  %%nuke_bands = find ( (fk_wide < 100) | (fk_wide >= 300 ) );
  %%nuke = [ nuke_lines(:); nuke_bands(:) ];
  nuke = [ nuke_lines(:) ];
  indsNuke = unique ( nuke(:) );
  indsNuke = indsNuke ( (indsNuke >= 1) & (indsNuke <= length(indsWide)) );

  %% ----- replace all data <100Hz, and >300Hz with Gaussian noise ----------
  DebugPrintf ( 1, "\n----- Identified lines in %s: -----\n", uvar.tsIn.IFO );
  DebugPrintf ( 1, "Line-frequencies: %f Hz\n", fk_wide ( inds_lines ) );

  if ( uvar.plotSpectrum )
    figure(); clf; hold on;
    plot ( fk_wide, abs ( xk_wide ) / sqrt(T), "+-", "color", "blue" ); legend ( uvar.tsIn.IFO );
    plot ( fk_wide, sqrt(Sn_wide), "o;sqrt(SX);", "color", "green" );

    NindsNuke = length ( indsNuke );
    draws = normrnd ( 0, 1, 2, NindsNuke );
    %%Sn0 = (8.2e-24)^2;
    Sn0 = Sn_wide ( indsNuke );
    noise = sqrt(T/4 * Sn0 ) .* ( draws(1,:) + I * draws(2,:) );
    plot ( fk_wide ( indsNuke ), abs ( xk_wide ( indsNuke )) / sqrt(T), "o", "color", "red" );
    xk_wide ( indsNuke ) = noise;
    plot ( fk_wide ( indsNuke ), abs ( xk_wide ( indsNuke )) / sqrt(T), "x", "color", "red" );

    %%draws = normrnd ( 0, 1, 2, NindsNuke );
    grid on;
    hold off;
  endif

  %% ----- whitened and overwhitened sFTs (with nuked lines) ----------
  xkW_wide  = xk_wide ./ sqrt(Sn_wide);
  xkOW_wide = xk_wide ./ Sn_wide;

  %% ----- turn SFTs back into timeseries ----------
  ts   = freqBand2TS ( fk_wide, xk_wide,   uvar.fMin, uvar.fMax, fSamp ); ts.ti   += tStart;
  tsW  = freqBand2TS ( fk_wide, xkW_wide,  uvar.fMin, uvar.fMax, fSamp ); tsW.ti  += tStart;
  tsOW = freqBand2TS ( fk_wide, xkOW_wide, uvar.fMin, uvar.fMax, fSamp ); tsOW.ti += tStart;

  assert ( (ts.ti == tsW.ti) && (ts.ti == tsOW.ti) );

  tsOut.ti   = ts.ti;
  tsOut.xi   = ts.xi;
  tsOut.xiW  = tsW.xi;
  tsOut.xiOW = tsOW.xi;
  tsOut.IFO = uvar.tsIn.IFO;

  return;
endfunction
