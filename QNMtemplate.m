function tsQNM = QNMtemplate ( t0GPS, A, phi0, f0, tau, ts, symmetric = false )
  ## return QNM template time-series of params {t0GPS, A, phi0, f0, tau} on given time-series 'ts'

  tsQNM = ts;
  tsQNM.xi(:) = 0;

  Dt_i = ts.ti - (t0GPS - ts.epoch);
  indsRingdown = find ( Dt_i >= 0 );
  indsRingup = [];
  if ( symmetric )
    indsRingup = find ( Dt_i < 0 );
  endif
  if ( isempty ( indsRingdown ) && isempty ( indsRingup ) )
    return;
  endif

  if ( !isempty ( indsRingdown ) )
    Dt_pos = Dt_i ( indsRingdown );
    tsQNM.xi ( indsRingdown )  = A * e.^(- Dt_pos / tau ) .* cos ( 2*pi * f0 * Dt_pos + phi0 );
  endif
  if ( !isempty ( indsRingup ) )
    Dt_neg = Dt_i ( indsRingup );
    tsQNM.xi ( indsRingup )  = A * e.^( Dt_neg / tau ) .* cos ( 2*pi * f0 * Dt_neg + phi0 );
  endif

  return;

endfunction
