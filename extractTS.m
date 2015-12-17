#!/usr/bin/octave -q

function [ts, tsW, tsOW, psd] = extractTS ( varargin )
  global debugLevel = 1;

  uvar = parseOptions ( varargin,
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"Twindow", "real,strictpos,scalar", 10 },	%% time-window +- to extract around the event
                        {"lineSigma", "real,positive,scalar", 5},	%% sigma deviations to indentify 'lines' in spectrum
                        {"lineWidth", "real,positive,scalar", 0.1},	%% +- width in Hz to zero around 'lines'
                        {"plotResults", "bool", false },
                        {"simulate", "bool", false },
                        {"shiftL", "real,positive,scalar", 7.1e-3},	%% time-shift to apply to L1 data wrt to H1
                        {"RngMedWindow", "real,positive,scalar", 300 }  %% window size to use for rngmed-based PSD estimation
                      );
  assert ( uvar.fMax > uvar.fMin );

  resDir = "extractTS-Results";
  [status, msg] = mkdir ( resDir );
  assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );

  %% load frequency-domain data from SFTs:
  fNy = 1370;	%% 2x1370Hz sampling, enough to allow for resolved ~7.3ms time-shift ~ 20bins
  fnames = {"./Data/H-1_H1_1800SFT_ER8-1126257832-1800.sft"; "./Data/L-1_L1_1800SFT_ER8-1126258841-1800.sft" };
  IFO = {"H1"; "L1"};
  tOffs = { 0, uvar.shiftL };	%% delay {H1,L1} data by this amount, respectively
  scaleFact = { 1, -1 };	%% invert L1 data to be in phase with H1

  if ( uvar.simulate )
    extraLabel = "-sim";
  else
    extraLabel = "";
  endif
  bname = sprintf ( "freq%.0fHz-%.0fHz-GPS%.0fs+-%.0fs-ls%.1f-lw%.1f-rng%.0f-shiftL%.2fms%s",
                    uvar.fMin, uvar.fMax, uvar.tCenter, uvar.Twindow, uvar.lineSigma, uvar.lineWidth, uvar.RngMedWindow, uvar.shiftL * 1e3, extraLabel);

  for X = 1:length(fnames)
    bnameX = sprintf ( "%s-%s", IFO{X}, bname );
    psd_fname = sprintf ( "%s/PSD-%s.dat", resDir, bnameX );
    ts_fname = sprintf ( "%s/TS-%s.dat", resDir, bnameX );
    %% ---------- check if TS results for this parameters already exist: re-use if yes ----------
    if ( length ( glob ( { psd_fname; ts_fname } ) ) == 2 )
      DebugPrintf (2, "%s: Re-using previous TS results '%s'\n", funcName(), bnameX );
      dat = load ( psd_fname );
      psd{X}.fk = (dat(:,1))';
      psd{X}.Sn = (dat(:,2))';
      dat = load ( ts_fname );
      ts{X}.ti   = tsW{X}.ti = tsOW{X}.ti = (dat(:,1))';
      ts{X}.xi   = (dat(:,2))';
      tsW{X}.xi  = (dat(:,3))';
      tsOW{X}.xi = (dat(:,4))';
    else
      %% ---------- otherwise: Read SFT frequency-domain data ----------
      DebugPrintf (2, "%s: Extracting TS from SFTs\n", funcName() );
      sft = readSFT ( fnames{X} );
      t0 = sft.header.epoch.gpsSeconds;
      f0 = sft.header.f0;
      assert ( strcmp ( IFO{X}, sft.header.IFO ) );
      Tsft = sft.header.Tsft; df = 1/Tsft;
      Nfreq = length ( sft.SFTdata );
      f1 = f0 + (Nfreq-1) * df;
      fk0{X} = f0 : df : f1;
      if ( uvar.simulate )
        sqrtSn = 8e-24;
        sigma = sqrt(Tsft)/2 * sqrtSn;
        xk0{X} = normrnd ( 0, sigma, 1, Nfreq ) + I * normrnd ( 0, sigma, 1, Nfreq );
      else
        xk0{X} = sft.SFTdata(:,1) + I * sft.SFTdata(:,2);
      endif
      assert ( length(fk0) == length(xk0) );

      %% ---------- extract frequency band of interest [fMin,fMax] as a timeseries ----------
      tsBand0 = freqBand2TS ( fk0{X}, xk0{X}, uvar.fMin - 1, uvar.fMax + 1, fNy );	%% 1Hz extra-band for median PSD estimates later
      tsBand0.ti += t0 + tOffs{X};	%% label times by correct epoch, shift by detector-specific offset
      tsBand0.xi *= scaleFact{X};		%% apply detector-specific scale factor to time-series

      %% ---------- truncate timeseries to [ tCenter - dT, tCenter + dT ] ----------
      indsTrunc = find ( (tsBand0.ti >= (uvar.tCenter - uvar.Twindow)) & (tsBand0.ti <= (uvar.tCenter + uvar.Twindow )) );
      tsBand.ti = tsBand0.ti ( indsTrunc );
      tsBand.xi = tsBand0.xi ( indsTrunc );

      %% ---------- compute PSD on short timeseries, nuke lines, extract 'physical' frequency band, and whiten + overwhitened TS ----------
      [psd{X}, ts{X}, tsW{X}, tsOW{X}] = whitenTS ( tsBand, uvar.fMin, uvar.fMax, uvar.lineSigma, uvar.lineWidth, uvar.RngMedWindow );

      fid = fopen ( psd_fname, "wb" );
      fprintf ( fid, "%%%% %18s %16s\n", "freq [Hz]", "SX [1/Hz]" );
      fprintf ( fid, "%16.9f  %16.9g\n", [psd{X}.fk', psd{X}.Sn']' );
      fclose(fid);

      fid = fopen ( ts_fname, "wb" );
      fprintf ( fid, "%%%% %18s %16s %16s %16s\n", "ti [GPS s]", "xi", "xi/sqrtSX", "xi/SX" );
      fprintf ( fid, "%18.9f %16.9g %16.9g %16.9g\n", [ts{X}.ti', ts{X}.xi', tsW{X}.xi', tsOW{X}.xi']' );
      fclose(fid);
    endif %% if no previous results re-used


    ts{X}.IFO = tsW{X}.IFO = tsOW{X}.IFO = psd{X}.IFO = IFO{X};

    if ( uvar.plotResults )
      ft = FourierTransform ( ts{X}.ti, ts{X}.xi );
      figure(); clf;
      plot ( ft.fk, abs(ft.xk)/sqrt(uvar.Twindow), "-", psd{X}.fk, sqrt(psd{X}.Sn), "o" );
      xlim ( [uvar.fMin, uvar.fMax] );
      xlabel ("Freq [Hz]");
      ylabel ("sqrt(Sn)");
      title ( bnameX );
      fname = sprintf ( "%s/%s.tex", resDir, bnameX);
      ezprint ( fname, "width", 512 );
    endif

  endfor %% X

  return;
endfunction
