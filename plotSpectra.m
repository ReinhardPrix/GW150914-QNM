function [H_psd, H_spect] = plotSpectra ( multiTS, multiPSD, xrange = [0, 2e3] )

  psd0_H1 = load ( "./Data/lalinferencemcmc-0-H1L1-1126259462.39-0H1-PSD.dat" );
  psd0_L1 = load ( "./Data/lalinferencemcmc-0-H1L1-1126259462.39-0L1-PSD.dat");

  numIFOs = length ( multiTS );
  assert ( length(multiPSD) == numIFOs );

  H_psd = figure(); clf;
  for X = 1 : numIFOs
    ts = multiTS{X};
    psd = multiPSD{X};

    dt = 1 / ts.fSamp;
    T  = dt * length(ts.xi);
    IFO = ts.IFO;

    if ( strcmp ( IFO, "H1" ) )
      psd0 = psd0_H1;
    elseif ( strcmp ( IFO, "L1" ) )
      psd0 = psd0_L1;
    else
      psd0 = [];
    endif

    if ( !isfield ( ts, "fk" ) )
      binsFreqRange = binRange ( ts.fMin, ts.fMax, psd.fk );
      %% compute narrow-banded FFT of input TS
      win = tukeywin ( length(ts.xi), 0.1 )';
      ts.xi .*= win;
      xFFT = dt * fft ( ts.xi );
      ts.xk = xFFT ( binsFreqRange );
      ts.fk = psd.fk ( binsFreqRange );
      Sn = psd.Sn ( binsFreqRange );
      ts.xkW = ts.xk ./ sqrt ( Sn );
      ts.xkOW = ts.xk ./ Sn;
      multiTS{X} = ts;
    endif

    %% ----- plot PSD for detector X
    subplot ( numIFOs, 1, X, "align" ); hold on;
    semilogy ( ts.fk, abs ( ts.xk ) / sqrt(T), sprintf ( "+-;%s;", IFO), "color", "blue" );
    if ( !isempty( "psd0" ) )
      semilogy ( psd0(:,1), sqrt(psd0(:,2)), "x", "color", "magenta" );
    endif
    semilogy ( psd.fk, sqrt(psd.Sn), "o", "color", "green" );
    xlim ( xrange );
    ylim ( [ 1e-24, 1e-20 ] );
    grid on;
    legend ( "location", "NorthEast" );
    ylabel ( "sqrt(Sn) [Hz^(-1/2)]" );

  endfor %% X = 1:numIFOs
  xlabel ("Freq [Hz]");

  H_spect = figure(); clf;
  for X = 1 : numIFOs
    ts = multiTS{X};
    IFO = ts.IFO;

    %% ----- plot whitened and over-whitened spectra for detector X
    subplot ( 2, numIFOs, X, "align" );
    plot ( ts.fk, abs ( ts.xkW ), sprintf("+-;%s;",IFO), "color", "blue" );
    xlim ( xrange );
    ylim ( [ 0, 20 ] );
    grid on;
    legend ( "location", "NorthEast" );
    if ( X == 1 ) ylabel ( "|x|/sqrt(SX)" ); endif

    subplot ( 2, numIFOs, X + numIFOs );
    plot ( ts.fk, abs ( ts.xkOW ), sprintf("+-;%s;",IFO), "color", "blue" );
    xlim ( xrange );
    ylim ( [ 0, 2e24 ] );
    xlabel ("Freq [Hz]");
    grid on;
    legend ( "location", "NorthEast" );
    if ( X == 1 ) ylabel ( "|x|/SX" ); endif

  endfor %% X = 1 : numIFOs
  xlabel ("Freq [Hz]");

  return;

endfunction %% plotSpectra()
