function [TS, FT] = freqBand2TS ( ft, fMin, fMax, fSamp )

  %% handy shortcuts
  fk = ft.fk(:);
  xk = ft.xk(:);
  IFO = ft.IFO;
  epoch = ft.epoch;
  df = mean ( diff ( fk ) );

  %% extract frequency-band
  inds0 = binRange ( fMin, fMax, fk );

  %% window FT bins before inv-FFTing:
  win = tukeywin ( length ( inds0 ), 0.1 );
  xkWin = xk(inds0) .* win;

  %% ========== place this band into a full spectrum including negative frequencies ==========
  Nfreq = round ( fSamp / df );
  %% taken from 'fftshift()':
  fk1 = df * [ - ( ceil ( (Nfreq-1)/2 ) : (-1) : 1 ), 0, ( 1 : floor ( (Nfreq-1)/2) ) ];
  assert ( Nfreq == length(fk1) );
  %% ----- positive frequencies ----------
  xk1 = zeros ( size(fk1) );
  inds1P = binRange ( fMin, fMax, fk1 );
  assert ( length(inds0) == length(inds1P) );
  xk1(inds1P) = xkWin;
  %% ----- negative frequencies ----------
  inds1N = binRange ( -fMax, -fMin, fk1 );
  assert ( length(inds0) == length(inds1N) );
  xk1(inds1N) = conj ( flipdim ( xkWin ) );	%% mirror-image and complex-conjugate
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
  win = tukeywin ( Nsamp, 0.01 );
  xiWin = TS.xi(:) .* win(:);
  TS.xi = xiWin';

  %% ===== also return narrow-banded original Fourier spectrum in [fMin, fMax] =====
  FT.fk = fk ( inds0 );
  FT.xk = xk ( inds0 );
  FT.IFO = IFO;
  FT.epoch = epoch;

  return;

endfunction
