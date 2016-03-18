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

function [TS, FT] = freqBand2TS ( ft, fMin, fMax, fSamp )

  %% handy shortcuts
  fk = ft.fk(:);
  xk = ft.xk(:);
  IFO = ft.IFO;
  epoch = ft.epoch;
  df = mean ( diff ( fk ) );

  %% extract frequency-band
  bins0 = binRange ( fMin, fMax, fk );

  %% window FT bins before inv-FFTing:
  win = tukeywin ( length ( bins0 ), 0.1 );
  xkWin = xk(bins0) .* win;

  %% ========== place this band into a full spectrum including negative frequencies ==========
  Nfreq = round ( fSamp / df );
  %% taken from 'fftshift()':
  fk1 = df * [ - ( ceil ( (Nfreq-1)/2 ) : (-1) : 1 ), 0, ( 1 : floor ( (Nfreq-1)/2) ) ];
  assert ( Nfreq == length(fk1) );
  %% ----- positive frequencies ----------
  xk1 = zeros ( size(fk1) );
  bins1P = binRange ( fMin, fMax, fk1 );
  assert ( length(bins0) == length(bins1P) );
  xk1(bins1P) = xkWin;
  %% ----- negative frequencies ----------
  bins1N = binRange ( -fMax, -fMin, fk1 );
  assert ( length(bins0) == length(bins1N) );
  xk1(bins1N) = conj ( flipdim ( xkWin ) );	%% mirror-image and complex-conjugate
  %% ----- inv-FFT back into time-domain ----------
  TS.xi = Nfreq * df * ifft ( ifftshift ( xk1 ) );
  dt = 1 / fSamp;
  Nsamp = length(TS.xi);	%% == Nfreq
  TS.ti = dt * [ 0 : (Nsamp-1) ];
  %% check that we constructed a real-valued timeseries
  err = max ( abs(imag( TS.xi)) ./ abs(real(TS.xi) ) );
  assert ( err < 1e-6 );
  TS.xi = real ( TS.xi );
  TS.epoch = ft.epoch;
  TS.IFO = ft.IFO;

  %% window final TS to avoid switch on/off artifacts:
  win = tukeywin ( Nsamp, 0.1 );
  xiWin = TS.xi(:) .* win(:);
  TS.xi = xiWin';

  %% ===== also return narrow-banded original Fourier spectrum in [fMin, fMax] =====
  FT.fk = fk ( bins0 );
  FT.xk = xk ( bins0 );
  FT.IFO = IFO;
  FT.epoch = epoch;

  return;

endfunction
