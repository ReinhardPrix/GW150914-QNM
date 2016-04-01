## Copyright (C) 2016 Reinhard Prix
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

function tsOut = injectQNMs ( tsIn, injectionSources = [] )
  %% tsInj = injectQNMs ( tsIn, injectionSources )
  %% inject QNM signal(s) into given timeseries 'ts'

  global debugLevel = 1;

  numInjections = length ( injectionSources );
  tsOut = tsIn;

  for l = 1 : numInjections
    inj_l = injectionSources(l);

    skyCorr = inj_l.skyCorr;
    X1 = find ( strcmp ( tsIn.IFO, {skyCorr.IFO} ) ); assert ( !isempty(X1) && (length(X1) == 1) );
    timeShift = skyCorr(X1).timeShift;
    ampFact   = skyCorr(X1).ampFact;
    t0    = inj_l.t0 - timeShift;
    A_eff = ampFact * inj_l.A

    t0InjOffs = t0 - tsIn.epoch;
    iStart_l = min ( find ( (tsIn.ti - t0InjOffs) >= 0 ) );	%% find earliest QNM time within scope
    if ( isempty ( iStart_l ) )
      DebugPrintf ( 2, "Dropping injection %d in %s: t0Inj = %.6f s > tEnd = %.6f s\n", l, tsIn.IFO, tsIn.epoch + t0InjOffs, tsIn.epoch + max ( tsIn.ti ) );
    else
      DebugPrintf ( 2, "Starting injection %d in %s at t_i(%d) = %.6f s, for t0 = %.6f s: {A = %g, phi0 = %g, f0 = %g Hz, tau = %g ms ... ",
                    l, tsIn.IFO, iStart_l, tsIn.epoch + tsIn.ti(iStart_l), tsIn.epoch + t0InjOffs, inj_l.A, inj_l.phi0, inj_l.f0, inj_l.tau  * 1e3 );
      tsQNM_l  = QNMtemplate ( t0, A_eff, inj_l.phi0, inj_l.f0, inj_l.tau, tsIn, (smooth = true) );	%% include a 'smooth' ramp-up to avoid injecting discontinuities
      tsOut.xi += tsQNM_l.xi;
      DebugPrintf ( 2, "done.\n");
    endif
  endfor %% for l = 1 : numInjections

  return;

endfunction %% injectQNMs()
