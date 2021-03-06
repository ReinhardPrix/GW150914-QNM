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

function plotContours ( resV, resCommon, select = [], plotMarkers = [] )
  global debugLevel = 1;

  Nsearches = length ( resV );
  ff0 = resCommon.ff0;
  ttau = resCommon.ttau;

  maxlog10B = 0;
  for l = 1 : Nsearches
    log10B_l = log10 ( resV(l).BSG );
    if ( log10B_l > maxlog10B )
      maxlog10B = log10B_l;
    endif
  endfor

  figure (); clf;
  hold on;
  if ( !isempty(plotMarkers) )
    for i = 1 : length ( plotMarkers )
      leg = sprintf ( ";%s;", plotMarkers(i).name );
      plot ( plotMarkers(i).f0, plotMarkers(i).tau * 1e3, leg, "marker", "o", "markersize", 5 );
    endfor
  endif

  colors = {[0,128,0]./256, [128,0,128]./256, [220,20,60]./256, [255, 165, 0]./256 };
  counter = 0;
  for l = 1 : Nsearches
    if ( !isempty(select) && all ( l != select ) ) continue; endif
    counter ++;
    log10B_l = log10 ( resV(l).BSG );
    DebugPrintf ( 1, "t0 = tM + %.1fms: log10BSG = %.2g, SNR = %.2g\n", resV(l).tOffs * 1e3, log10B_l, resV(l).AmpMP.SNR );
    color_l = colors { mod (counter - 1, length(colors) ) + 1 };
    if ( !isempty ( resV(l).isoConf2 ) )
      [C, H] = contour ( ff0, ttau * 1e3, resV(l).BSG_f0_tau, resV(l).isoConf2 * [ 1, 1 ] );
    endif
    set ( H, "linecolor", color_l, "linewidth", 5, "linestyle", "--" );
    leg = sprintf ( "%.1fms", resV(l).tOffs * 1e3 );
    cl = clabel ( C, H, "FontSize", 12, "Color", color_l);
    set ( cl, "string", "" );
    %%set ( cl(1), "string", leg );
    set ( cl(end), "string", leg );
    plot ( resV(l).lambdaMP.f0, resV(l).lambdaMP.tau * 1e3, "*", 'markeredgecolor', color_l, "markersize", 6 );
  endfor
  xlabel ( "Freq [Hz]" );
  ylabel ( "tau [ms]" );
  ylim ( [ 0, 14.5 ] );
  hold off;
  grid on;

  return;

endfunction %% plotContours()
