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

function plotSummary ( resV, resCommon )

  %% ----- prepare quantities to be plotted versus QNM start-time
  tMerger = resCommon.tMerger;
  Nsteps = length ( resV );
  BSG = f0_MPE = f0_lerr = f0_uerr = taums_MPE = taums_lerr = taums_uerr = zeros ( 1, Nsteps );
  for i = 1 : Nsteps
    BSG(i)        = resV(i).BSG;

    f0_MPE(i)     = resV(i).f0_est.MPE;
    f0_lerr(i)    = resV(i).f0_est.lerr;
    f0_uerr(i)    = resV(i).f0_est.uerr;

    taums_MPE(i)  = resV(i).tau_est.MPE * 1e3;
    taums_lerr(i) = resV(i).tau_est.lerr * 1e3;
    taums_uerr(i) = resV(i).tau_est.uerr * 1e3;
  endfor

  tOffsVms = [resV.tOffs] * 1e3;
  tOffs_Range = [ (min ( tOffsVms(:) )), (max ( tOffsVms(:) )) + 0.1 ];
  taums_Range = [ (min ( resCommon.ttau(:) )), (max ( resCommon.ttau(:))) + 1e-4 ] * 1e3;
  f0_Range    = [ (min ( resCommon.ff0(:) )),  (max ( resCommon.ff0(:) )) + 0.1 ];

  %% ===== plot summary 1 ==========
  figure (); clf;

  %% ----- plot log10BSG(tOffs)
  subplot ( 2, 2, 1, "align" );
  xrange = tOffs_Range;
  yrange = [ -1, 20 ];
  plot ( tOffsVms, log10([resV.BSG]), "-o" );
  xlim ( xrange );
  ylim ( yrange );
  line ( xrange, 0, "linestyle", "-", "linewidth", 3 );
  line ( xrange, 1, "linestyle", ":", "linewidth", 3 );
  line ( 0 * [1,1], yrange, "linestyle", "-", "linewidth", 2 );
  xlim ( xrange );
  ylabel ("log10<BSG>");
  grid on;

  %% ----- plot f0_MPE(tOffs)
  subplot ( 2, 2, 2, "align" );
  errorbar ( tOffsVms, f0_MPE, f0_lerr, f0_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( f0_Range );
  ylabel ("f0 [Hz]");

  %% ----- plot SNR(tOffs)
  subplot ( 2, 2, 3, "align" ); hold on;
  plot ( tOffsVms, [[resV.AmpMP].SNR], "-o;MP;" );
  plot ( tOffsVms, [[resV.AmpML].SNR], "-x;ML;" );
  grid on;
  xlim ( xrange );
  ylabel ("SNR");
  xlabel ( sprintf ( "%.6f s + tOffs [ms]", tMerger) );

  %% ----- plot tau_MPE(tOffs)
  subplot ( 2, 2, 4, "align" );
  errorbar ( tOffsVms, taums_MPE, taums_lerr, taums_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( taums_Range );
  ylabel ("tau [ms]");
  xlabel ( sprintf ( "%.6f s + tOffs [ms]", tMerger) );

  fname = sprintf ( "%s-summary.pdf", resCommon.bname );
  ezprint ( fname, "width", 512 );


  return;

endfunction
