#!/usr/bin/octave -q

function [ts, ft, psd] = extractTSfromSFT ( varargin )
  global debugLevel = 1;
  global use_v2 = false;

  uvar = parseOptions ( varargin,
                        {"SFTpath", "char,vector" },
                        {"fMin", "real,strictpos,scalar", 100 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"Twindow", "real,strictpos,scalar", 10 },	%% time-window +- to extract around the event
                        {"lineSigma", "real,positive,scalar", 5},	%% sigma deviations to indentify 'lines' in spectrum
                        {"lineWidth", "real,positive,scalar", 0.1},	%% +- width in Hz to zero around 'lines'
                        {"plotSpectrum", "bool", false },
                        {"RngMedWindow", "real,positive,scalar", 300 },  %% window size to use for rngmed-based PSD estimation
                        {"fSamp", "real,positive,scalar", 2000*2 },	%% sampling rate of output timeseries
                        {"useBuffer", "bool", true}			%% re-use timeseries-data files if found
                      );
  assert ( uvar.fMax > uvar.fMin );

  pieces = strsplit ( uvar.SFTpath, "/" );
  sftfname = pieces { end };
  sftbname = strrep ( sftfname, ".sft", "");
  bname = sprintf ( "%s-freq%.0fHz-%.0fHz-fSamp%.0fHz-GPS%.0fs+-%.0fs-ls%.1f-lw%.1f-rng%.0f",
                    sftbname, uvar.fMin, uvar.fMax, uvar.fSamp, uvar.tCenter, uvar.Twindow, uvar.lineSigma, uvar.lineWidth, uvar.RngMedWindow );

  resDir = "extractTS-Results";
  [status, msg] = mkdir ( resDir );
  assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );

  psd_fname = sprintf ( "%s/%s.psd", resDir, bname );
  ts_fname  = sprintf ( "%s/%s.ts",  resDir, bname );
  ft_fname  = sprintf ( "%s/%s.ft",  resDir, bname );
  %% extract IFO name from SFT name: only works for SFT-name compliant SFTs
  pieces = strsplit ( sftfname, {"-", "_"} );
  IFO = pieces{3};
  assert ( (length ( IFO ) == 2) && isalpha(IFO(1)) && isdigit(IFO(2)) );

  %% ---------- check if TS results for this parameters already exist: re-use if yes ----------
  if ( uvar.useBuffer && (length ( glob ( { psd_fname; ts_fname; ft_fname } ) ) == 3) )
    DebugPrintf (2, "%s: Re-using previous TS results '%s'\n", funcName(), bname );

    dat = load ( psd_fname );
    psd.fk = (dat(:,1))';
    psd.Sn = (dat(:,2))';
    psd.IFO = IFO;

    dat = load ( ft_fname );
    ft.fk   = (dat(:,1))';
    ft.xk   = (dat(:,2) + I * dat(:,3))';
    ft.xkW  = (dat(:,4) + I * dat(:,5))';
    ft.xkOW = (dat(:,6) + I * dat(:,6))';
    ft.IFO  = IFO;

    dat = load ( ts_fname );
    ts.ti   = (dat(:,1))';
    ts.xi   = (dat(:,2))';
    ts.xiW  = (dat(:,3))';
    ts.xiOW = (dat(:,4))';
    ts.IFO  = IFO;

    epoch = unique ( dat(:,5) );
    ts.epoch = epoch;
    psd.epoch = epoch;
    ft.epoch = epoch;
  else
    %% ---------- otherwise: Read SFT frequency-domain data ----------
    DebugPrintf (2, "%s: Extracting TS from SFT '%s'\n", funcName(), uvar.SFTpath );
    sft = readSFT ( uvar.SFTpath );
    epoch = sft.header.epoch.gpsSeconds;
    f0 = sft.header.f0;
    assert ( strcmp ( IFO, sft.header.IFO ) );
    Tsft = sft.header.Tsft; df = 1/Tsft;
    Nfreq = length ( sft.SFTdata );
    f1 = f0 + (Nfreq-1) * df;
    fk0 = f0 : df : f1;
    xk0 = sft.SFTdata(:,1) + I * sft.SFTdata(:,2);
    assert ( length(fk0) == length(xk0) );
    ft0.fk = fk0;
    ft0.xk = xk0;
    ft0.IFO = IFO;
    ft0.epoch = epoch;

    %% ---------- extract frequency band of interest [fMin,fMax] as a timeseries ----------
    tsBand0 = freqBand2TS ( ft0, uvar.fMin, uvar.fMax, uvar.fSamp );

    %% ---------- truncate timeseries to [ tCenter - dT, tCenter + dT ] ----------
    indsTrunc = find ( (tsBand0.ti >= (uvar.tCenter - tsBand0.epoch - uvar.Twindow)) & (tsBand0.ti <= (uvar.tCenter - tsBand0.epoch + uvar.Twindow )) );
    tsBand.ti = tsBand0.ti ( indsTrunc );
    tsBand.xi = tsBand0.xi ( indsTrunc );
    tsBand.IFO = IFO;
    tsBand.epoch = tsBand0.epoch;

    %% ---------- compute PSD on short timeseries, nuke lines, extract 'physical' frequency band, and whiten + overwhitened TS ----------
    if ( use_v2 )
      [ts, ft, psd] = whitenTS_v2 ( "tsIn", tsBand, ...
                                    "fMin", uvar.fMin, "fMax", uvar.fMax, ...
                                    "lineSigma", uvar.lineSigma, "lineWidth", uvar.lineWidth, ...
                                    "plotSpectrum", uvar.plotSpectrum );
    else
      [ts, ft, psd] = whitenTS ( "tsIn", tsBand,
                                 "fMin", uvar.fMin, "fMax", uvar.fMax,
                                 "lineSigma", uvar.lineSigma, "lineWidth", uvar.lineWidth,
                                 "RngMedWindow", uvar.RngMedWindow,
                                 "plotSpectrum", uvar.plotSpectrum );
    endif
    if ( uvar.plotSpectrum )
      title ( bname );
      fname = sprintf ( "%s/%s-spectrum.pdf", resDir, bname);
      ezprint ( fname, "width", 512 );
    endif

    %% ---------- store results for potential future re-use ----------
    fid = fopen ( psd_fname, "wb" ); assert ( fid != -1, "Failed to open '%s' for writing\n", psd_fname );
    fprintf ( fid, "%%%%%14s %16s\n", "freq [Hz]", "SX [1/Hz]" );
    fprintf ( fid, "%16.9f  %16.9g\n", [psd.fk', psd.Sn']' );
    fclose(fid);

    fid = fopen ( ts_fname, "wb" ); assert ( fid != -1, "Failed to open '%s' for writing\n", ts_fname );
    fprintf ( fid, "%%%%%14s %16s %16s %16s %16s\n", "ti [offs s]", "xi", "xi/sqrtSX", "xi/SX", "epoch [GPS s]" );
    epochV = ts.epoch * ones ( size ( ts.ti ) );	%% stupid way, but ensure we keep epoch separate for 'offsets' ti to avoid roundoff problems
    fprintf ( fid, "%16.9f %16.9g %16.9g %16.9g %16.9f\n", [ts.ti', ts.xi', ts.xiW', ts.xiOW', epochV']' );
    fclose(fid);

    fid = fopen ( ft_fname, "wb" ); assert ( fid != -1, "Failed to open '%s' for writing\n", ft_fname );
    fprintf ( fid, "%%%%%14s %16s %16s %16s %16s %16s %16s\n", "fk [Hz]", "Re(xk)", "Im(xk)", "Re(xk/sqrtSX)", "Im(xk/sqrtSX)", "Re(xk/SX)", "Im(xk/SX)" );
    fprintf ( fid, "%16.9f %16.9g %16.9g %16.9g %16.9g %16.9g %16.9g\n",
              [ ft.fk', real(ft.xk'), imag(ft.xk'), real(ft.xkW'), imag(ft.xkW'), real(ft.xkOW'), imag(ft.xkOW') ]' );
    fclose(fid);

  endif %% if no previous results re-used

  return;
endfunction
