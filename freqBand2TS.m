function TS = freqBand2TS ( fk, xk, fMin, fMax, fSamp )

  inds0 = find ( (fk >= fMin) & (fk <= fMax) );

  %% place this band into a full spectrum including negative frequencies
  df = mean ( diff ( fk ) );

  %% fk1 = -fNyquist : df : fNyquist;
  Nfreq = round ( fSamp / df );
  %% taken from 'fftshift()':
  fk1 = [ -(ceil((Nfreq-1)/2):-1:1)*df, 0, (1:floor((Nfreq-1)/2))*df ];
  xk1 = zeros ( size(fk1) );
  inds1P = find ( (fk1 >= fMin) & (fk1 <= fMax) );
  assert ( length(inds0) == length(inds1P) );
  xk1(inds1P) = xk(inds0);

  inds1N = find ( (fk1 >= -fMax) & (fk1 <= -fMin) );
  assert ( length(inds0) == length(inds1N) );
  xk1(inds1N) = conj ( xk( flipdim (inds0) ) );	%% mirror-image and complex-conjugate

  TS = FourierTransformInv ( fk1, xk1 );

  %% check that we constructed a real-valued timeseries
  err = max ( abs(imag( TS.xi)) ./ abs(real(TS.xi) ) );
  assert ( err < 1e-6 );
  TS.xi = real ( TS.xi );

  return;

endfunction
