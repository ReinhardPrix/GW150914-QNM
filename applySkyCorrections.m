function multiTS = applySkyCorrections ( multiTS, skyCorr )
  %% multiTS = applySkyCorrections ( multiTS, skyCorr )
  %% function to apply sky-dependent time-shift and antenna-pattern amplitude factors
  %% 'skyCorr' to the given multi-IFO time-series

  numIFOs = length ( multiTS );
  assert ( length(skyCorr) == numIFOs );

  dt = 1 / multiTS{1}.fSamp;

  for X = 1 : numIFOs
    IFO = multiTS{X}.IFO;
    X1 = find ( strcmp ( IFO, {skyCorr.IFO} ) ); assert ( !isempty(X1) && (length(X1) == 1) );
    shiftBins = round ( skyCorr(X1).timeShift / dt );
    ampFact   = skyCorr(X1).ampFact;
    shift_eff = shiftBins * dt; assert ( abs(shift_eff - skyCorr(X1).timeShift) < 1e-6 );	%% check we're within 0.001ms

    %% there's 2 ways to do the time-shifts: add shift-value to ti, or shift bins in xi
    %% for now we use the second method, as we're going to sum the time-shifted ts{X} together and
    %% compute a single Bayes-factor.
    multiTS{X}.xi   = ampFact * circshift ( multiTS{X}.xi,   [0, shiftBins] );
    if ( isfield ( multiTS{X}, "xiW" ) )
      multiTS{X}.xiW  = ampFact * circshift ( multiTS{X}.xiW,  [0, shiftBins] );
    endif
    if ( isfield ( multiTS{X}, "xiOW" ) )
      multiTS{X}.xiOW = ampFact * circshift ( multiTS{X}.xiOW, [0, shiftBins] );
    endif

  endfor %% X = 1 : numIFOs

  return;

endfunction %% applySkyCorrections()
