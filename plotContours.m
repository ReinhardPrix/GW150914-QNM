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

function plotContours ( in, select = [] )
  global iFig0 = 0;
  global tMergerOffs
  global tEvent;
  global f0GR;
  global taumsGR;

  Nsteps = length ( in );
  if ( isempty ( select ) )
    select = [ 1 : Nsteps ];
  endif

  maxlog10B = 0;
  for i = select
    log10Bi = log10 ( in{i}.BSG );
    if ( log10Bi > maxlog10B )
      maxlog10B = log10Bi;
    endif
  endfor

  figure ( iFig0 + 5 ); clf;
  hold on;
  he = errorbar ( f0GR.val, taumsGR.val, f0GR.lerr, f0GR.uerr, taumsGR.lerr, taumsGR.uerr, "~>.r;IMR;" );
  for i = select
    log10Bi = log10 ( in{i}.BSG );
    color_i = "black"; %% (1 - log10Bi / maxlog10B) * [ 1, 1, 1];
    if ( !isempty ( in{i}.isoConf2 ) )
      [C, H] = contour ( in{i}.ff0, in{i}.ttau * 1e3, in{i}.posterior2D, in{i}.isoConf2 * [ 1, 1 ] );
    endif
    set ( H, "linecolor", color_i, "linewidth", 2 );
    leg = sprintf ( "tM + %.1fms", (in{i}.tGPS - tEvent - tMergerOffs) * 1e3 );
    cl = clabel ( C, H, "FontSize", 12, "Color", color_i);
    set ( cl, "string", "" );
    %%set ( cl(1), "string", leg );
    set ( cl(end), "string", leg );
    %%plot ( in{i}.f0_MPE2, in{i}.tau_MPE2 * 1e3, "x;MPE;", 'markeredgecolor', color_i, "markersize", 5 );
  endfor
  xlabel ( "Freq [Hz]" );
  ylabel ( "tau [ms]" );
  hold off;
  grid on;

  fname = sprintf ( "%s-contours.pdf", in{1}.bname );
  ezprint ( fname, "width", 512 );

  return;

endfunction
