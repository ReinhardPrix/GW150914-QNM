%% estimate single-sided noise PSD 'psd(f)', and return whitened and over-whitened TS
%% including an automated line-nuking algorithm: 'lines' are identified as > lineSigma deviations in power
function [ tsOut, ftOut, psd ] = whitenTS_v2 ( varargin )
  global debugLevel = 1;
  global iFig0 = 0;

  uvar = parseOptions ( varargin,
                        {"ftIn", "struct" },
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"Twindow", "real,strictpos,scalar", 10 },	%% time-window +- to extract around the event
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"fSamp", "real,positive,scalar", 2000*2 },	%% sampling rate of output timeseries
                        {"plotSpectrum", "bool", false }
                      );

  %% handy shortcuts
  IFO = uvar.ftIn.IFO;

  %% turn SFT into a full-band time-series at output sampling rate
  fMinSFT = min ( uvar.ftIn.fk );
  fMaxSFT = max ( uvar.ftIn.fk );
  ts0 = freqBand2TS ( uvar.ftIn, fMinSFT, fMaxSFT, uvar.fSamp );

  %% identify time indices of truncated time-segement [tCenter - Twindow, tCenter + Twindow]
  samplesTrunc = find ( (ts0.ti >= (uvar.tCenter - uvar.ftIn.epoch - uvar.Twindow)) & (ts0.ti <= (uvar.tCenter - uvar.ftIn.epoch + uvar.Twindow )) );
  tsTrunc = ts0;
  tsTrunc.ti = ts0.ti ( samplesTrunc );
  tsTrunc.xi = ts0.xi ( samplesTrunc );
  %% re-adjust epoch, make 'offsets' ti start from 0:
  offs0 = tsTrunc.ti(1);
  tsTrunc.epoch += offs0;
  tsTrunc.ti    -= offs0;

  Nsamples = length ( samplesTrunc );
  dt = mean ( diff ( tsTrunc.ti ) );
  T = Nsamples * dt;

  %% ---------- estimate PSD on full-length (1800s) timeseries using pwelch() using segments of length 'samplesTrunc' ----------
  window = length ( samplesTrunc );
  overlap = 0.5;
  NsFFT = length ( window );
  [Sn, fkPSD] = pwelch ( ts0.xi, window, overlap, NsFFT, uvar.fSamp, 'onesided' );	%% computes single-sided PSD
  %% limit to original SFT frequency range
  binsFullPSD  = binRange ( fMinSFT, fMaxSFT, fkPSD );
  psd.fk    = fkPSD ( binsFullPSD )';
  psd.Sn    = Sn (binsFullPSD )';
  psd.IFO   = IFO;
  psd.epoch = ts0.epoch;

  %% ---------- compute full spectrum of truncated time-series for whitening ----------
  ftTrunc = FourierTransform ( tsTrunc.ti, tsTrunc.xi );
  ftTrunc.IFO = IFO;
  ftTrunc.epoch = tsTrunc.epoch;
  %% limit to original SFT frequency range
  binsFull  = binRange ( fMinSFT, fMaxSFT, ftTrunc.fk );
  ftTrunc.fk = ftTrunc.fk ( binsFull );
  ftTrunc.xk = ftTrunc.xk ( binsFull );
  %% PSD bins should agree exactly with data spectrum bins
  assert ( length ( binsFullPSD ) == length ( binsFull ) );
  err = max ( abs ( ftTrunc.fk(:) - psd.fk(:) ));
  assert ( err < 1e-6 );

  %% ----- whitened and overwhitened full-frequency range spectra ----------
  ftTrunc.xkW  = ftTrunc.xk ./ sqrt ( psd.Sn );
  ftTrunc.xkOW = ftTrunc.xk ./ psd.Sn;
  ftTrunc.xkU  = ftTrunc.xkW / sqrt(T/2);	%% normalized by rms=> unit variance

  %% ----- turn spectra back into narrow-banded timeseries ----------
  ft0    = ftTrunc;
  ts     = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, uvar.fSamp );

  ft0.xk = ftTrunc.xkW;
  tsW    = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, uvar.fSamp );

  ft0.xk = ftTrunc.xkOW;
  tsOW   = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, uvar.fSamp );
  assert ( (ts.ti == tsW.ti) && (ts.ti == tsOW.ti) );

  tsOut      = ts;
  tsOut.xiW  = tsW.xi;
  tsOut.xiOW = tsOW.xi;

  %% ----- return + plot narrow-banded spectra ----------
  ftOut      = ftTrunc;
  bins_out   = binRange ( uvar.fMin, uvar.fMax, ftTrunc.fk );
  ftOut.fk   = ftTrunc.fk ( bins_out );
  ftOut.xk   = ftTrunc.xk ( bins_out );
  ftOut.xkW  = ftTrunc.xkW ( bins_out );
  ftOut.xkU  = ftTrunc.xkU ( bins_out );	%% normalized by rms=> unit variance
  ftOut.xkOW = ftTrunc.xkOW ( bins_out );

  %% ---------- plot spectra over narrow-banded output frequency range ----------
  if ( uvar.plotSpectrum )
    if ( strcmp ( IFO, "H1" ) )
      psd0 = load ( "./Data/lalinferencemcmc-0-H1L1-1126259462.39-0H1-PSD.dat" );
    elseif ( strcmp ( IFO, "L1" ) )
      psd0 = load ( "./Data/lalinferencemcmc-0-H1L1-1126259462.39-0L1-PSD.dat");
    endif
    iIFO = ifelse ( strcmp ( IFO, "H1" ), 0, 1 );
    figure ( 5 * iFig0 + 1 + iIFO ); clf;

    subplot ( 3, 1, 1 ); hold on;
    plot ( ftOut.fk, abs ( ftOut.xk ) / sqrt(T), "+-", "color", "blue" ); legend ( IFO );
    if ( exist ( "psd0" ) )
      plot ( psd0(:,1), sqrt(psd0(:,2)), "x;lalinference;", "color", "magenta" );
    endif
    plot ( psd.fk, sqrt(psd.Sn), "o;sqrt(SX);", "color", "green" );
    xlim ( [uvar.fMin, uvar.fMax] );
    ylim ( [0, 5e-23] );
    grid on;
    hold off;

    subplot ( 3, 1, 2 );
    plot ( ftOut.fk, abs ( ftOut.xkU ), "+-", "color", "blue" ); legend ( "xk/rms" );
    xlim ( [uvar.fMin, uvar.fMax] );
    ylim ( [ 0, 5 ] );

    subplot ( 3, 1, 3 );
    plot ( ftOut.fk, abs ( ftOut.xkOW ), "+-", "color", "blue" ); legend ( "xk/Sn" );
    xlim ( [uvar.fMin, uvar.fMax] );
  endif

  return;

endfunction %% whitenTS_v2()
