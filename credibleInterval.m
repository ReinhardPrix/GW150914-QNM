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
function x_est = credibleInterval ( x, posterior_x, confidence = 0.9 )

  assert ( isvector ( x ) && isvector(posterior_x) && (length ( x ) == length ( posterior_x )) );

  normx = sum ( posterior_x(:) );
  posterior1D = posterior_x / normx;

  try
    [pIso, delta, INFO, OUTPUT] = fzero ( @(pIso)  sum ( posterior1D ( find ( posterior1D >= pIso ) ) ) - confidence, ...
                                          [ min(posterior1D), max(posterior1D) ], ...
                                          optimset ( "TolX", 1e-4 )
                                        );
    assert ( INFO == 1 );
  catch
    DebugPrintf (0, "fzero() failed\n");
    x_est = [];
    return;
  end_try_catch

  x_MP = x ( find ( posterior1D == max(posterior1D) ) );

  inds0 = find ( posterior1D >= pIso );
  i_min = min(inds0);
  i_max = max(inds0);
  %% check if these are respective closest to p_iso
  if ( (i_min > 1) && abs(posterior1D(i_min-1) - pIso) < abs(posterior1D(i_min) - pIso) )
    i_min --;
  endif
  if ( (i_max < length(posterior1D)) && abs(posterior1D(i_max+1) - pIso) < abs(posterior1D(i_max) - pIso) )
    i_max ++;
  endif

  x_est = struct ( "MPE", x_MP, "lerr", x_MP - x(i_min), "uerr", x(i_max) - x_MP, "pIso", pIso * normx );
  return;
endfunction
