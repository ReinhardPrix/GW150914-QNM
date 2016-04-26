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

function multiTS = injectQNMs ( multiTS, injectionSources = [] )
  %% multiTS = injectQNMs ( multiTS, injectionSources )
  %% inject QNM signal(s) 'injectionSources' into multi-IFO timeseries 'multiTS'

  global debugLevel = 1;

  numInjections = length ( injectionSources );
  numIFOs = length ( multiTS );

  for X = 1 : numIFOs

    tsX = multiTS{X};

    for l = 1 : numInjections
      inj_l = injectionSources(l);
      skyCorr = inj_l.skyCorr;
      X1 = find ( strcmp ( tsX.IFO, {skyCorr.IFO} ) ); assert ( !isempty(X1) && (length(X1) == 1) );
      timeShift = skyCorr(X1).timeShift;
      ampFact   = skyCorr(X1).ampFact;

      t0Eff    = inj_l.t0 - timeShift;	%% time-shifted start-time
      AEff     = ampFact * inj_l.A;	%% antenna-pattern corrected amplitude

      t0Eff_offs = t0Eff - tsX.epoch;
      iStart_l = min ( find ( (tsX.ti - t0Eff_offs) >= 0 ) );	%% find earliest QNM time within scope
      if ( isempty ( iStart_l ) )
        DebugPrintf ( 2, "Dropping injection %d in %s: t0Inj = %.6f s > tEnd = %.6f s\n", l, tsX.IFO, tsX.epoch + t0Eff_offs, tsX.epoch + max ( tsX.ti ) );
      else
        DebugPrintf ( 2, "Starting injection %d in %s at t_i(%d) = %.6f s, for t0 = %.6f s: {A = %g, phi0 = %g, f0 = %g Hz, tau = %g ms ... ",
                      l, tsX.IFO, iStart_l, tsX.epoch + tsX.ti(iStart_l), tsX.epoch + t0Eff_offs, inj_l.A, inj_l.phi0, inj_l.f0, inj_l.tau  * 1e3 );
        tsQNM_l  = QNMtemplate ( t0Eff, AEff, inj_l.phi0, inj_l.f0, inj_l.tau, tsX, (ringup = true) );	%% include a 'ringup' to avoid injecting discontinuities
        tsX.xi  += tsQNM_l.xi;
        DebugPrintf ( 2, "done.\n");
      endif
    endfor %% for l = 1 : numInjections

    multiTS{X} = tsX;

  endfor %% X = 1 : numIFOs

  return;

endfunction %% injectQNMs()
