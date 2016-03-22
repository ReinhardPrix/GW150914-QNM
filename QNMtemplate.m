function tsQNM = QNMtemplate ( t0GPS, A, phi0, f0, tau, ts )
  ## return QNM template time-series of params {t0GPS, A, phi0, f0, tau} on given time-series 'ts'
  xi = [];

  Dt_i = ts.ti - (t0GPS - ts.epoch);
  indsRingdown = find ( Dt_i >= 0 );
  if ( isempty ( indsRingdown ) )
    return;
  endif

  Dt_pos = Dt_i ( indsRingdown );
  tsQNM.xi = A * e.^(- Dt_pos / tau ) .* cos ( 2*pi * f0 * Dt_pos + phi0 );
  tsQNM.epoch = t0GPS + min(Dt_pos);
  tsQNM.ti = Dt_pos - min(Dt_pos);	%% timesteps start from 0

  return;

endfunction
