## Copyright (C) 2015 Reinhard Prix
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

function plotSnapshot ( in, ts, tRange = [0.41, 0.455], plotMarkers = [] )
  global iFig0 = 0;
  global tMergerOffs
  global tEvent;

  f0 = unique ( in.ff0 );
  f0Range = [ min(f0), max(f0) ];

  tau = unique ( in.ttau );
  tauRange = [ min(tau), max(tau) ];

  f0_est  = in.f0_est;
  tau_est = in.tau_est;

  figure ( iFig0 + 3 ); clf;

  %% ----- posterior (f0, tau)
  subplot ( 2, 2, 1 );
  hold on;
  colormap ("jet");
  surf ( in.ff0, in.ttau * 1e3, in.posterior2D ); view(2); shading("interp"); %% colorbar("location", "NorthOutside");
  plot3 ( in.f0_MPE2, in.tau_MPE2 * 1e3,  1.1 * in.posteriorMax, "marker", "x", "markersize", 2, "linewidth", 2, "color", "white" );
  xlim ( f0Range );
  ylim ( tauRange * 1e3 );
  %%xlabel ("f0 [Hz]");
  ylabel ("tau [ms]");

  if ( !isempty(plotMarkers) )
    for i = 1 : length ( plotMarkers )
      plot3 ( plotMarkers(i).f0, plotMarkers(i).tau * 1e3, 1.1 * in.posteriorMax, "marker", "o", "markersize", 3, "color", "white" );
    endfor
  endif

  if ( !isempty(in.isoConf2))
    [C, H] = contour ( in.ff0, in.ttau * 1e3, in.posterior2D, in.isoConf2 * [ 1, 1 ] );
    set ( H, "linecolor", "white", "linewidth", 1 );
  endif
  hold off;

  %% ----- posterior(f0)
  subplot ( 2, 2, 3 );
  plot ( f0, in.posterior_f0, "linewidth", 2 );
  grid on;
  yrange = ylim();
  if ( !isempty ( in.f0_est ) )
    line ( [in.f0_est.MPE, in.f0_est.MPE], yrange );
    line ( [in.f0_est.MPE - in.f0_est.lerr, in.f0_est.MPE + in.f0_est.uerr], [in.f0_est.pIso, in.f0_est.pIso] );
  endif
  xlim ( f0Range );
  ylim ( yrange );
  xlabel ("f0 [Hz]");
  %%ylabel ("pdf(f0)");
  set ( gca(), "yticklabel", {} );

  %% ----- posterior (tau)
  subplot ( 2, 2, 2 );
  plot ( in.posterior_tau, tau * 1e3, "linewidth", 2 );
  xrange = xlim();
  if ( !isempty ( in.tau_est ) )
    line ( xrange, [in.tau_est.MPE, in.tau_est.MPE]*1e3 );
    line ( [in.tau_est.pIso, in.tau_est.pIso], [in.tau_est.MPE - in.tau_est.lerr, in.tau_est.MPE + in.tau_est.uerr]* 1e3 );
  endif
  xlim ( xrange );
  ylim ( tauRange * 1e3 );
  grid on;
  %%label ("pdf(tau)");
  ylabel ("tau [ms]");
  set ( gca(), "xticklabel", {} );

  %% ----- timeseries {h(t), template}
  subplot ( 2, 2, 4 );
  hold on;
  colors = { "red", "blue" };
  %% plot data timeseries:
  for X = 1 : length ( ts )
    sleg = sprintf (";OW[%s] ;", ts{X}.IFO );
    %% prepare per-detector timeseries for matching: adapt L1 data to be phase-coherent with H1
    shiftBins = 0;
    if ( strcmp ( ts{X}.IFO, "L1" ) )
      shiftL1 = 7.0e-3;
      dt = mean( diff ( ts{X}.ti ) );
      shiftBins = round ( shiftL1 / dt );
      shiftL1_eff = shiftBins * dt;
      assert ( abs(shiftL1_eff - shiftL1) < 1e-6 );
      %% adjust time-offset and relative (antenna-pattern) amplitudes
      ts{X}.ti += shiftL1_eff;
      ts{X}.xi   *= -1;
      ts{X}.xiW  *= -1;
      ts{X}.xiOW *= -1;
    endif

    plot ( ts{X}.epoch + ts{X}.ti - tEvent, ts{X}.xiOW * ts{X}.SX_GR, sleg, "linewidth", 2, "color", colors{X} );
  endfor

  tOffs = in.tGPS - tEvent;
  tOffsFromM = tOffs - tMergerOffs;

  %% MPE template
  tmpl_MPE = QNMtemplate ( in.tGPS, in.A_MPE2, in.phi0_MPE2, in.f0_MPE2, in.tau_MPE2, ts{1} );
  plot ( tmpl_MPE.ti - tEvent, tmpl_MPE.xi, ";MPE ;", "linewidth", 4, "color", "black" );
  legend ( "location", "NorthEast");
  yrange = [-1.7e-21, 1.5e-21 ];
  line ( (in.tGPS - tEvent) * [ 1, 1 ],   yrange, "linestyle", "--", "linewidth", 2 );
  line ( [ tMergerOffs, tMergerOffs], yrange, "linestyle", "-", "linewidth", 2 );
  if ( !isempty ( tRange ) )
    xlim ( tRange );
  else
    xlim ( [tOffs - 0.01, tOffs + 5 * in.tau_MPE2 ] );
  endif
  ylim ( yrange );
  xlabel ( sprintf ( "%.0f + tOffs [s]", tEvent) );
  text ( min(xlim()) - 0.2 * abs(diff(xlim())), 0, "h(t)" );

  textOffs = sprintf ( "tOffs = %.5fs = tM + %.2fms", tOffs, tOffsFromM * 1e3 );
  title ( textOffs );
  textB = sprintf ( "log10(BSG) = %.2g\nSNR0 = %.2g", log10 ( in.BSG ), in.SNR_MPE2 );
  x0 = tOffs + 0.02 * abs(diff(xlim()));
  y0 = min(yrange) + 0.2*abs(diff(yrange));
  text ( x0, y0, textB );
  ylim ( yrange );
  grid on;
  hold off;

  fname = sprintf ( "%s.pdf", in.bname );
  ezprint ( fname, "width", 512 );

endfunction %% plotSnapshot()
