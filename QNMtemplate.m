function tsQNM = QNMtemplate ( t0GPS, A, phi0, f0, tau, ts, ringup = false )
  ## return QNM template time-series of params {t0GPS, A, phi0, f0, tau} on given time-series 'ts'

  tsQNM = ts;
  tsQNM.xi(:) = 0;

  Dt_i = ts.ti - (t0GPS - ts.epoch);
  indsRingdown = find ( Dt_i >= 0 );
  indsRingup = [];
  if ( ringup )
    indsRingup = find ( Dt_i < 0 );
  endif
  if ( isempty ( indsRingdown ) && isempty ( indsRingup ) )
    return;
  endif

  if ( !isempty ( indsRingdown ) )
    Dt_pos = Dt_i ( indsRingdown );
    tsQNM.xi ( indsRingdown )  = A * e.^(- Dt_pos / tau ) .* cos ( 2*pi * f0 * Dt_pos + phi0 );
  endif
  %% if 'ringup' option given: include a 10x faster 'ringup' to avoid discontinuous start-up
  if ( !isempty ( indsRingup ) )
    Dt_neg = Dt_i ( indsRingup );
    tsQNM.xi ( indsRingup )  = A * e.^( Dt_neg / (tau/10) ) .* cos ( 2*pi * f0 * Dt_neg + phi0 );
  endif

  return;

endfunction

%!demo
%! fSamp = 4e3;
%! dt = 1/fSamp;
%! ts.ti = [0 : 100] * dt;
%! ts.epoch = 0;
%! t0 = 20*dt; A = 1; phi0 = pi/4; f0 = 251; tau = 4e-3;
%! tsQNM = QNMtemplate ( t0, A, phi0, f0, tau, ts, (ringup = false) );
%! tsQNMru = QNMtemplate ( t0, A, phi0, f0, tau, ts, (ringup = true) );
%! figure(); clf; hold on;
%! plot ( tsQNM.ti, tsQNM.xi, "o;QNM;", tsQNMru.ti, tsQNMru.xi, "-x;ringup+QNM;" );
%! line ( t0 * [1,1], ylim, "linestyle", "--", "linewidth", 3 );
%! grid on;
%! hold off;
%! title ( sprintf ( "t0 = %.3f s; A = %g, phi0 = %.1f, f0 = %.1f Hz, tau = %.1f ms", t0, A, phi0, f0, tau * 1e3 ) );
%! xlabel ( "t [s]" );
%! ezprint ( "QNMtemplate.pdf", "width", 256 );
