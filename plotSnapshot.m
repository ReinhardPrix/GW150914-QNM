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

function plotSnapshot ( in, tRange = [0.41, 0.455], plotMarkers = [] )
  global iFig0 = 0;

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
  markers = {".", "."};
  %% plot data timeseries:
  for X = 1 : length ( in.ts )
    sleg = sprintf (";OW[%s] ;", in.ts{X}.IFO );
    Dti_ms = 1e3 * ( in.ts{X}.epoch + in.ts{X}.ti - in.tMerger );
    plot ( Dti_ms, in.ts{X}.xiOW * in.ts{X}.SX_GR, sleg, "linewidth", 2, "color", colors{X}, "marker", markers{X} );
  endfor

  tOffs_ms = in.tOffs * 1e3;
  tau_ms = in.tau_MPE2 * 1e3;
  %% MPE template in H1
  tmpl_MPE = QNMtemplate ( in.t0GPS, in.A_MPE2, in.phi0_MPE2, in.f0_MPE2, in.tau_MPE2, in.ts{1} );	%% refer to IFO 1, assumed H1
  Dti_ms = (tmpl_MPE.epoch + tmpl_MPE.ti - in.tMerger) * 1e3;
  plot( Dti_ms, tmpl_MPE.xi, ";MPE ;", "linewidth", 4, "color", "black" );
  legend ( "location", "NorthEast");
  yrange = [-1.7e-21, 1.5e-21 ];
  xrange = [tOffs_ms - 10, tOffs_ms + 5 * tau_ms ];
  line ( [ 0, 0 ],   yrange, "linestyle", "--", "linewidth", 2 );
  line ( tOffs_ms * [1,1], yrange, "linestyle", "-", "linewidth", 2 );
  xlabel ( sprintf ( "%.6f s + tOffs [ms]", in.tMerger ) );
  %%text ( min(xlim()) - 0.2 * abs(diff(xlim())), 0, "h(t)" );
  ylabel ( "h(t)");

  textB = sprintf ( "log10(BSG) = %.2g\nSNR0 = %.2g\nt0 = +%.1f ms", log10 ( in.BSG ), in.SNR_MPE2, tOffs_ms );
  x0 = tOffs_ms + 0.02 * abs(diff(xrange));
  y0 = min(yrange) + 0.2*abs(diff(yrange));
  text ( x0, y0, textB );
  x0 = 0 + 0.02 * diff(xrange);
  y0 = min(yrange) + 0.07*abs(diff(yrange));
  text ( x0, y0, "tMerger" );

  xlim ( xrange );
  ylim ( yrange );
  grid on;
  hold off;

  fname = sprintf ( "%s.pdf", in.bname );
  ezprint ( fname, "width", 512 );

endfunction %% plotSnapshot()
