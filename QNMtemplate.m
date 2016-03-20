function tsQNM = QNMtemplate ( t0, A, phi0, f0, tau, ts )
  ## return QNM template time-series of params {t0, A, phi0, f0, tau} on given time-series 'ts'
  xi = [];

  Dt_i = ts.ti - (t0 - ts.epoch);
  indsRingdown = find ( Dt_i >= 0 );
  if ( isempty ( indsRingdown ) )
    return;
  endif

  Dt_pos = Dt_i ( indsRingdown );
  xi = A * e.^(- Dt_pos / tau ) .* cos ( 2*pi * f0 * Dt_pos + phi0 );
  ti = t0 + Dt_pos;

  tsQNM.ti = ti;
  tsQNM.xi = xi;
  return;

endfunction
