function psd = estimatePSD ( ts, sampleRange )
  %% psd = estimatePSD ( ts, sampleRange )
  %% estimate PSD on full-length timeseries *excluding* the analysis segment 'sampleRange' = [i_t0, i_t1]
  %%
  %% Note: the 'sampleRange' is used as the welch window-length, so the returned frequency
  %% bins are df = 1/duration(sampleRange)

  assert ( isvector ( sampleRange ) && (length(sampleRange) == 2) );

  samplesSeg = [ sampleRange(1) : sampleRange(2) ];
  windowLen = length ( samplesSeg );
  overlap = 0.5;
  NsFFT = windowLen;
  xiBefore = ts.xi ( 1 : (sampleRange(1)-1) );
  xiAfter  = ts.xi ( sampleRange(2) + 1 : end );
  [SnBefore, fkBefore] = pwelch ( xiBefore, windowLen, overlap, NsFFT, ts.fSamp, 'onesided' );	%% computes single-sided PSD
  [SnAfter,  fkAfter]  = pwelch ( xiAfter,  windowLen, overlap, NsFFT, ts.fSamp, 'onesided' );	%% computes single-sided PSD
  assert ( (length(fkBefore) == length(fkAfter)) && all ( (fkBefore == fkAfter)(:) ) );
  fkSn = fkBefore';
  Sn   = 0.5 * ( SnBefore + SnAfter )';

  psd.IFO = ts.IFO;
  psd.epoch = ts.epoch;
  psd.fk = fkSn;
  psd.Sn = Sn;

  return;

endfunction %% estimatePSD()
