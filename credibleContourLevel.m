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

function zConf = credibleContourLevel ( posterior2D, confidence = 0.9 )
  %% return 2D contour iso-pdf level containing 'confidence' probability of the total

  assert  ( ismatrix ( posterior2D ) );
  assert ( (confidence > 0) && (confidence <1 ) );

  normC = sum ( posterior2D (:) );
  posterior2D ./= normC;
  try
    [zConf0, delta, INFO, OUTPUT] = fzero ( @(zIso)  sum ( posterior2D ( find ( posterior2D >= zIso ) )(:) ) - confidence, ...
                                            [ min(posterior2D(:)), max(posterior2D(:)) ], ...
                                            optimset ( "TolX", 1e-4 )
                                          );
    assert ( INFO == 1 );
    zConf = [];
  catch
    printf ("fzero() failed\n");
    zConf = [];
    return;
  end_try_catch

  zConf = normC * zConf0;
  return;

endfunction
