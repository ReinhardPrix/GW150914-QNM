function [multiTS, multiPSD] = whitenNarrowBandTS ( multiTS, multiPSD, FreqRange )
  %% multiTS = whitenNarrowBandTS ( multiTS, multiPSD, FreqRange )

  numIFOs = length ( multiTS );
  assert ( length ( multiPSD ) == numIFOs );
  assert ( isvector ( FreqRange ) && (length(FreqRange) == 2) );
  fMin = min(FreqRange); fMax = max(FreqRange);

  for X = 1 : numIFOs

    ts = multiTS{X};
    psd = multiPSD{X};

    fSamp = ts.fSamp;
    dt = 1 / fSamp;

    %% ----- extract frequency band 'FreqRange' ----------
    binsFreqRange = binRange ( fMin, fMax, psd.fk );

    %% truncate psd to frequency band
    psd.fk = psd.fk ( binsFreqRange );
    psd.Sn = psd.Sn ( binsFreqRange );

    %% compute narrow-banded FFT of input TS
    win = tukeywin ( length(ts.xi), 0.1 )';
    ts.xi .*= win;
    xFFT = dt * fft ( ts.xi );
    xk = xFFT ( binsFreqRange );

    ft.IFO   = ts.IFO;
    ft.epoch = ts.epoch;
    ft.fk    = psd.fk;
    ft.xk    = xk;

    %% whiten and over-whiten data
    ftW = ft;
    ftW.xk   = xk ./ sqrt ( psd.Sn );

    ftOW = ft;
    ftOW.xk  = xk ./ psd.Sn;

    %% ----- turn spectra back into narrow-banded timeseries ----------
    ts   = freqBand2TS ( ft,   fMin, fMax, fSamp );
    tsW  = freqBand2TS ( ftW,  fMin, fMax, fSamp );
    tsOW = freqBand2TS ( ftOW, fMin, fMax, fSamp );
    assert ( all((ts.ti == tsW.ti)(:)) && all((ts.ti == tsOW.ti)(:)) );

    %% ----- store return values
    multiPSD{X} = psd;

    multiTS{X} = ts;
    multiTS{X}.xiW  = tsW.xi;
    multiTS{X}.xiOW = tsOW.xi;

    %% also store spectra, for plotting
    multiTS{X}.fk   = ft.fk;
    multiTS{X}.xk   = ft.xk;
    multiTS{X}.xkW  = ftW.xk;
    multiTS{X}.xkOW = ftOW.xk;

  endfor %% X = 1 : numIFOs

  return;

endfunction %% whitenNarrowBandTS()

