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

%% return MP estimate 'x_est' with fields 'MPE', 'lower', 'upper', and 'pIso'
%% ie: maximum-posterior estimate x_est.MPE
%% 'confidence'-credible interval [x_est.lower, x_est.upper], and
%% corresponding iso-posterior value x_est.pIso
function x_est = credibleInterval ( posterior, confidence = 0.9 )

  assert ( isstruct ( posterior ) && isvector([posterior.x]) && isvector ([posterior.px]) && (length([posterior.x]) == length([posterior.px])) );
  assert ( sum ( posterior.px(:) ), 1, 1e-6 );

  x      = posterior.x;
  prob_x = posterior.px;
  try
    [pIso, delta, INFO, OUTPUT] = fzero ( @(pIso)  sum ( prob_x ( find ( prob_x >= pIso ) ) ) - confidence, ...
                                          [ min(prob_x), max(prob_x) ]
                                        );
    assert ( INFO == 1 );
  catch
    DebugPrintf (0, "fzero() failed\n");
    x_est = [];
    return;
  end_try_catch

  [val, iMaxP ] = max ( prob_x(:) );
  x_MP = [posterior.x] ( iMaxP );

  inds0 = find ( prob_x >= pIso );
  i_min = min(inds0);
  i_max = max(inds0);
  %% check if these are respective closest to p_iso
  if ( (i_min > 1) && abs(prob_x(i_min-1) - pIso) < abs(prob_x(i_min) - pIso) )
    i_min --;
  endif
  if ( (i_max < length(prob_x)) && abs(prob_x(i_max+1) - pIso) < abs(prob_x(i_max) - pIso) )
    i_max ++;
  endif

  x_est = struct ( "MPE", x_MP, "lerr", x_MP - x(i_min), "uerr", x(i_max) - x_MP, "pIso", pIso );
  return;
endfunction
