%% estimate single-sided noise PSD 'psd(f)', and return whitened and over-whitened TS
%% including an automated line-nuking algorithm: 'lines' are identified as > lineSigma deviations in power
function [ tsOut, ftOut, psd ] = whitenTS_v2 ( varargin )
  global debugLevel = 1;

  uvar = parseOptions ( varargin,
                        {"tsIn", "struct" },
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"lineSigma", "real,positive,scalar", 5},	%% sigma deviations to indentify 'lines' in spectrum
                        {"lineWidth", "real,positive,scalar", 0.1},	%% +- width in Hz to zero around 'lines'
                        {"fSamp", "real,positive,scalar", 2000*2 },	%% sampling rate of output timeseries
                        {"plotSpectrum", "bool", false }
                      );

  %% handy shortcuts
  epoch = uvar.tsIn.epoch;
  ti = uvar.tsIn.ti;
  xi = uvar.tsIn.xi;
  IFO = uvar.tsIn.IFO;

  tStart = min ( ti );
  ti -= tStart;
  epoch += tStart;
  dt = mean ( diff ( ti ) );
  fSamp = round ( 1/dt );
  T = max(ti) + dt;
  ft = FourierTransform ( ti, xi );
  ft.IFO = IFO;
  ft.epoch = epoch;

  ft = FourierTransform ( ti, xi );
  indsUse  = binRange ( uvar.fMin, uvar.fMax, ft.fk );
  fk_use   = ft.fk ( indsUse );
  xk_use = ft.xk ( indsUse );

  %% ---------- estimate PSD using pwelch() over time-segemnents ----------
  Nsamples = length ( ti );
  window = 2^(ceil( log2 ( sqrt(Nsamples) ) ) );	%% 'default' window size
  overlap = 0.5;
  NsFFT = length ( ti ); %% use oversampling to get identical PSD freq-bins with input data-FT

  [Sn, fkPSD] = pwelch ( xi, window, overlap, NsFFT, fSamp, 'onesided' );	%% computes single-sided PSD
  indsUsePSD = binRange ( uvar.fMin, uvar.fMax, fkPSD );
  assert ( length(indsUse) == length(indsUsePSD) );
  psd.Sn = Sn ( indsUsePSD )';
  psd.fk = fkPSD ( indsUsePSD )';
  psd.IFO   = IFO;
  psd.epoch = epoch;

  assert ( fk_use(1), psd.fk(1), 1e-6 );
  assert ( fk_use(end), psd.fk(end), 1e-6 );

  if ( uvar.plotSpectrum )
    iFig = 4 + ifelse ( strcmp ( IFO, "H1" ), 1, 2 );
    figure(iFig); clf; hold on;
    plot ( fk_use, abs ( xk_use ) / sqrt(T), "+-", "color", "blue" ); legend ( IFO );
    plot ( psd.fk, sqrt(psd.Sn), "o;sqrt(SX);", "color", "green" );
    grid on;
    hold off;
    xlim ( [uvar.fMin, uvar.fMax] );
    ylim ( [0, 5e-23] );
  endif

  %% ----- whitened and overwhitened sFTs (with nuked lines) ----------
  xkW_use  = xk_use ./ sqrt(psd.Sn);
  xkOW_use = xk_use ./ psd.Sn;

  %% ----- turn SFTs back into timeseries ----------
  ftOut.fk = fk_use;
  ftOut.IFO = IFO;
  ftOut.epoch = epoch;
  ftOut.xk   = xk_use;
  ftOut.xkW  = xkW_use;
  ftOut.xkOW = xkOW_use;

  ft0 = ftOut;
  ts   = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, fSamp );
  ft0.xk = ftOut.xkW;
  tsW  = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, fSamp );
  ft0.xk = xkOW_use;
  tsOW = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, fSamp );

  assert ( (ts.ti == tsW.ti) && (ts.ti == tsOW.ti) );

  tsOut.ti   = ts.ti;
  tsOut.xi   = ts.xi;
  tsOut.xiW  = tsW.xi;
  tsOut.xiOW = tsOW.xi;
  tsOut.IFO =  IFO;
  tsOut.epoch = epoch;

  return;
endfunction
