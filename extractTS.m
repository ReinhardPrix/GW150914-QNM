#!/usr/bin/octave -q

function [ts, psd] = extractTS ( varargin )
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
                        {"RngMedWindow", "real,positive,scalar", 300 },  %% window size to use for rngmed-based PSD estimation
                        {"fSamp", "real,positive,scalar", 2000*2 }	%% sampling rate of output timeseries
                      );
  assert ( uvar.fMax > uvar.fMin );

  resDir = "extractTS-Results";
  [status, msg] = mkdir ( resDir );
  assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );

  %% load frequency-domain data from SFTs:
  fnames = {"./Data/H-1_H1_1800SFT_ER8-1126257832-1800.sft"; "./Data/L-1_L1_1800SFT_ER8-1126258841-1800.sft" };
  IFO = {"H1"; "L1"};

  if ( uvar.simulate )
    extraLabel = "-sim";
  else
    extraLabel = "";
  endif
  bname = sprintf ( "freq%.0fHz-%.0fHz-fSamp%.0fHz-GPS%.0fs+-%.0fs-ls%.1f-lw%.1f-rng%.0f%s",
                    uvar.fMin, uvar.fMax, uvar.fSamp, uvar.tCenter, uvar.Twindow, uvar.lineSigma, uvar.lineWidth, uvar.RngMedWindow, extraLabel);

  sideband = uvar.RngMedWindow * ( 1 / (2*uvar.Twindow)); 		%% extra frequency side-band for median PSD estimates later

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
      psd{X}.IFO = IFO{X};
      dat = load ( ts_fname );
      ts{X}.ti   = (dat(:,1))';
      ts{X}.xi   = (dat(:,2))';
      ts{X}.xiW  = (dat(:,3))';
      ts{X}.xiOW = (dat(:,4))';
      ts{X}.IFO  = IFO{X};
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
      tsBand0 = freqBand2TS ( fk0{X}, xk0{X}, uvar.fMin - sideband, uvar.fMax + sideband, uvar.fSamp );

      %% ---------- truncate timeseries to [ tCenter - dT, tCenter + dT ] ----------
      indsTrunc = find ( (tsBand0.ti >= (uvar.tCenter - t0 - uvar.Twindow)) & (tsBand0.ti <= (uvar.tCenter - t0 + uvar.Twindow )) );
      tsBand.ti = tsBand0.ti ( indsTrunc );
      tsBand.xi = tsBand0.xi ( indsTrunc );
      tsBand.IFO = IFO{X};

      %% ---------- compute PSD on short timeseries, nuke lines, extract 'physical' frequency band, and whiten + overwhitened TS ----------
      [psd{X}, ts{X}] = whitenTS ( tsBand, uvar.fMin, uvar.fMax, uvar.lineSigma, uvar.lineWidth, uvar.RngMedWindow );
      ts{X}.ti += t0;		%% label times by correct epoch

      fid = fopen ( psd_fname, "wb" );
      fprintf ( fid, "%%%% %18s %16s\n", "freq [Hz]", "SX [1/Hz]" );
      fprintf ( fid, "%16.9f  %16.9g\n", [psd{X}.fk', psd{X}.Sn']' );
      fclose(fid);

      fid = fopen ( ts_fname, "wb" );
      fprintf ( fid, "%%%% %18s %16s %16s %16s\n", "ti [GPS s]", "xi", "xi/sqrtSX", "xi/SX" );
      fprintf ( fid, "%18.9f %16.9g %16.9g %16.9g\n", [ts{X}.ti', ts{X}.xi', ts{X}.xiW', ts{X}.xiOW']' );
      fclose(fid);
    endif %% if no previous results re-used

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
