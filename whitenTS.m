%% estimate single-sided noise PSD 'psd(f)', and return whitened and over-whitened TS
%% including an automated line-nuking algorithm: 'lines' are identified as > lineSigma deviations in power
function [psd, ts, tsW, tsOW] = whitenTS ( tsIn, fMin, fMax, lineSigma = 5, lineWidth=0.1, RngMedWindow = 300 )
  fudge = 10 * eps;

  ti = tsIn.ti;
  tStart = min ( ti );
  dt = mean ( diff ( ti ) );
  T = max(ti) - tStart + dt;

  ft = FourierTransform ( tsIn.ti, tsIn.xi );
  fNyquist = max ( ft.fk );
  indsWide = find ( (ft.fk >=  fMin - 1) & (ft.fk <=  fMax + 1) );

  indsUse  = find ( (ft.fk(indsWide) >= fMin * (1 - fudge)) & (ft.fk(indsWide) <=  fMax * ( 1 + fudge ) ) );

  %% ---------- estimate PSD using running-median over frequency ----------
  %% Wiener-Khintchine theorm
  %% Sn_double = 1/T * mean ( periodo );
  window = 300;
  rngmedbias = 0.69817; %%
  lal;
  rngmedbias = XLALRngMedBias ( window );

  periodo = abs ( ft.xk(indsWide) ).^2;
  Sn_wide = 2/T * rngmed ( periodo, window ) / rngmedbias;	%% single-sided PSD
  psd.fk = ft.fk ( indsWide ( indsUse ) );
  psd.Sn = Sn_wide ( indsUse );

  %% ---------- automatically identify and nuke 'lines' ----------
  fk_wide   = ft.fk ( indsWide );
  xk_wide   = ft.xk ( indsWide );

  yk = xk_wide  ./ sqrt( T * Sn_wide / 2 );   %% normalized SFT:
  Pk_re = abs ( real ( yk ) );
  Pk_im = abs ( imag ( yk ) );
  inds_lines = find ( (Pk_re > lineSigma) | (Pk_im > lineSigma) );
  inukeWidth = round ( lineWidth * T );
  for i = 1: length(inds_lines)
    inds_nuke = [ inds_lines(i) - inukeWidth : inds_lines(i) + inukeWidth ];
    xk_wide ( inds_nuke ) = 0;
  endfor

  %% ----- whitened and overwhitened sFTs (with nuked lines) ----------
  xkW_wide  = xk_wide ./ sqrt(Sn_wide);
  xkOW_wide = xk_wide ./ Sn_wide;

  %% ----- turn SFTs back into timeseries ----------
  ts   = freqBand2TS ( fk_wide, xk_wide,   fMin, fMax, fNyquist ); ts.ti   += tStart;
  tsW  = freqBand2TS ( fk_wide, xkW_wide,  fMin, fMax, fNyquist ); tsW.ti  += tStart;
  tsOW = freqBand2TS ( fk_wide, xkOW_wide, fMin, fMax, fNyquist ); tsOW.ti += tStart;

  assert ( (ts.ti == tsW.ti) && (ts.ti == tsOW.ti) );
  return;
endfunction
