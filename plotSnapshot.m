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

function plotSnapshot ( res, resCommon, plotMarkers = [] )

  ff0  = resCommon.ff0;
  ttau = resCommon.ttau;

  ts = applySkyCorrections ( resCommon.multiTS, resCommon.skyCorr );

  tMerger = res.tMerger;

  f0 = unique ( ff0(:) );
  f0Range = [ min(f0), max(f0) ];

  tau = unique ( ttau(:) );
  tauRange = [ min(tau), max(tau) ];

  figure (); clf;

  %% ----- posterior (f0, tau)
  subplot ( 2, 2, 1 );
  hold on;
  colormap ("jet");
  surf ( ff0, ttau * 1e3, res.BSG_f0_tau ); view(2); shading("interp"); %% colorbar("location", "NorthOutside");
  plot3 ( res.lambdaMP.f0, res.lambdaMP.tau * 1e3,  1.1 * res.lambdaMP.BSG_MP, "marker", "x", "markersize", 2, "linewidth", 2, "color", "white" );
  xlim ( f0Range );
  ylim ( tauRange * 1e3 );
  %%xlabel ("f0 [Hz]");
  ylabel ("tau [ms]");

  if ( !isempty(plotMarkers) )
    for i = 1 : length ( plotMarkers )
      plot3 ( plotMarkers(i).f0, plotMarkers(i).tau * 1e3, 1.1 * res.lambdaMP.BSG_MP, "marker", "o", "markersize", 3, "linewidth", 2, "color", "white" );
    endfor
  endif

  if ( isfield ( res, "isoConf2") && !isempty(res.isoConf2))
    [C, H] = contour ( ff0, ttau * 1e3, res.BSG_f0_tau, res.isoConf2 * [ 1, 1 ] );
    set ( H, "linecolor", "white", "linewidth", 1 );
  endif
  hold off;

  %% ----- posterior(f0)
  subplot ( 2, 2, 3 );
  plot ( f0, res.posterior_f0.px, "linewidth", 2 );
  grid on;
  yrange = ylim();
  if ( isfield ( res, "f0_est" ) && !isempty ( res.f0_est ) )
    line ( [res.f0_est.MPE, res.f0_est.MPE], yrange );
    line ( [res.f0_est.MPE - res.f0_est.lerr, res.f0_est.MPE + res.f0_est.uerr], [res.f0_est.pIso, res.f0_est.pIso] );
  endif
  xlim ( f0Range );
  ylim ( yrange );
  xlabel ("f0 [Hz]");
  %%ylabel ("pdf(f0)");
  set ( gca(), "yticklabel", {} );

  %% ----- posterior (tau)
  subplot ( 2, 2, 2 );
  plot ( res.posterior_tau.px, tau * 1e3, "linewidth", 2 );
  xrange = xlim();
  if ( isfield ( res, "tau_est" ) && !isempty ( res.tau_est ) )
    line ( xrange, [res.tau_est.MPE, res.tau_est.MPE]*1e3 );
    if ( !isna ( res.tau_est.pIso ) )
      line ( [res.tau_est.pIso, res.tau_est.pIso], [res.tau_est.MPE - res.tau_est.lerr, res.tau_est.MPE + res.tau_est.uerr]* 1e3 );
    endif
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
  markers = {".", "o"};
  markersize = [ 5, 2 ];
  %% plot data timeseries:
  for X = 1 : length ( ts )
    sleg = sprintf (";OW[%s] ;", ts{X}.IFO );
    Dti_ms = 1e3 * ( ts{X}.epoch + ts{X}.ti - tMerger );
    plot ( Dti_ms, ts{X}.xiOW * resCommon.SXmid(X), sleg, "linewidth", 2, "color", colors{X}, "marker", markers{X}, "markersize", markersize(X) );
  endfor

  tOffs_ms = (res.t0 - res.tMerger ) * 1e3;
  tau_ms = res.lambdaMP.tau * 1e3;
  %% MPE template in H1
  tmpl_MPE = QNMtemplate ( res.t0, res.AmpMP.A, res.AmpMP.phi0, res.lambdaMP.f0, res.lambdaMP.tau, ts{1} );	%% refer to IFO 1, assumed H1
  Dti_ms = (tmpl_MPE.epoch + tmpl_MPE.ti - tMerger) * 1e3;
  plot( Dti_ms, tmpl_MPE.xi, ";MPE ;", "linewidth", 4, "color", "black" );
  legend ( "location", "NorthEast");
  yrange = [-1.8e-21, 1.7e-21 ];
  taumsMax = 1e3 * max ( tauRange );
  xrange = [tOffs_ms - taumsMax/2, tOffs_ms + 1.3 * taumsMax ];
  line ( [ 0, 0 ],   yrange, "linestyle", "-", "linewidth", 2 );
  line ( tOffs_ms * [1,1], yrange, "linestyle", "--", "linewidth", 3 );
  xlabel ( sprintf ( "H1:%.6f s + tOffs [ms]", tMerger ) );
  %%text ( min(xlim()) - 0.2 * abs(diff(xlim())), 0, "h(t)" );
  ylabel ( "h(t)");

  textB = sprintf ( "log10(BSG) = %.2g\nSNR0 = %.2g\ntM + %.1f ms", log10 ( res.BSG ), res.AmpMP.SNR, tOffs_ms );
  x0 = tOffs_ms + 0.02 * abs(diff(xrange));
  y0 = min(yrange) + 0.2*abs(diff(yrange));
  text ( x0, y0, textB );

  xlim ( xrange );
  ylim ( yrange );
  grid on;
  hold off;

endfunction %% plotSnapshot()
