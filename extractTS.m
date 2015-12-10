#!/usr/bin/octave -q

function [ts, tsW, tsOW, psd] = extractTS ( varargin )

  uvar = parseOptions ( varargin,
                        {"fMin", "real,strictpos,scalar", 200 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"tCenter", "real,strictpos,scalar", 1126259462 },
                        {"Twindow", "real,strictpos,scalar", 10 },	%% time-window +- to extract around the event
                        {"lineSigma", "real,positive,scalar", 5},	%% sigma deviations to indentify 'lines' in spectrum
                        {"lineWidth", "real,positive,scalar", 0.1},	%% +- width in Hz to zero around 'lines'
                        {"showPlots", "bool", false },
                        {"simulate", "bool", false },
                        {"RngMedWindow", "real,positive,scalar", 100 }	%% window size to use for rngmed-based PSD estimation
                      );
  assert ( uvar.fMax > uvar.fMin );

  %% load frequency-domain data from SFTs:
  fNy = 1370;	%% 2x1370Hz sampling, enough to allow for resolved ~7.3ms time-shift ~ 20bins
  fnames = {"H-1_H1_1800SFT_ER8-1126257832-1800.sft"; "L-1_L1_1800SFT_ER8-1126258841-1800.sft" };
  tOffs = { 0, +7.3e-3 };	%% delay L1 data by 7.3ms
  scaleFact = { 1, -1 };	%% invert L1 data to be in phase with H1

  if ( uvar.simulate )
    extraLabel = "-sim";
  else
    extraLabel = "";
  endif
  bname = sprintf ( "freq%.0fHz-%.0fHz-GPS%.0fs+-%.0fs%s", uvar.fMin, uvar.fMax, uvar.tCenter, uvar.Twindow, extraLabel);

  for X = 1:length(fnames)

    %% ---------- Read SFT frequency-domain data ----------
    sft = readSFT ( fnames{X} );
    t0 = sft.header.epoch.gpsSeconds;
    f0 = sft.header.f0;
    IFO{X} = sft.header.IFO;
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
    bnameX = sprintf ( "%s-%s", IFO{X}, bname );

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

    if ( uvar.showPlots )
      ft = FourierTransform ( ts{X}.ti, ts{X}.xi );
      figure(); clf;
      plot ( ft.fk, abs(ft.xk)/sqrt(uvar.Twindow), "-", psd{X}.fk, sqrt(psd{X}.Sn), "o" );
      xlim ( [uvar.fMin, uvar.fMax] );
      xlabel ("Freq [Hz]");
      ylabel ("sqrt(Sn)");
      title ( bnameX );
    endif

    psd_fname = sprintf ( "PSD-%s.dat", bnameX );
    fid = fopen ( psd_fname, "wb" );
    fprintf ( fid, "%.9f  %g\n", [psd{X}.fk', psd{X}.Sn']' );
    fclose(fid);

    ts_fname = sprintf ( "TS-%s.dat", bnameX );
    fid = fopen ( ts_fname, "wb" );
    fprintf ( fid, "%.9f  %g\n", [ts{X}.ti', ts{X}.xi']' );
    fclose(fid);

    tsW_fname = sprintf ( "TSW-%s.dat", bnameX );
    fid = fopen ( tsW_fname, "wb" );
    fprintf ( fid, "%.9f  %g\n", [tsW{X}.ti', tsW{X}.xi']' );
    fclose(fid);

    tsOW_fname = sprintf ( "TSOW-%s.dat", bnameX );
    fid = fopen ( tsOW_fname, "wb" );
    fprintf ( fid, "%.9f  %g\n", [tsOW{X}.ti', tsOW{X}.xi']' );
    fclose(fid);

  endfor %% X

  if ( uvar.showPlots )
    sleg1 = sprintf (";%s;", IFO{1} );
    sleg2 = sprintf (";%s;", IFO{2} );

    figure(); clf; hold on;
    plot ( ts{1}.ti - uvar.tCenter, ts{1}.xi, sleg1, "linewidth", 3, ts{2}.ti - uvar.tCenter, ts{2}.xi, sleg2, "linewidth", 3 );
    xlim ( [0.38, 0.46 ] );
    xlabel ("tGPS - tCenter [s]");
    ylabel ("Strain");
    title ( sprintf ( "TS - %s", bname ) );
    grid on;
    yrange = [-1.2, 1.2] * 1e-21;
    ylim ( yrange );
    hold off;

    figure(); clf; hold on;
    plot ( tsW{1}.ti - uvar.tCenter, tsW{1}.xi, sleg1, "linewidth", 3, tsW{2}.ti - uvar.tCenter, tsW{2}.xi, sleg2, "linewidth", 3 );
    xlim ( [0.38, 0.46 ] );
    xlabel ("tGPS - tCenter [s]");
    ylabel ("Strain/sqrt(SX(f))");
    title ( sprintf ( "TS-W-%s", bname ) );
    grid on;
    hold off;

    figure(); clf; hold on;
    plot ( tsOW{1}.ti - uvar.tCenter, tsOW{1}.xi, sleg1, "linewidth", 3, tsOW{2}.ti - uvar.tCenter, tsOW{2}.xi, sleg2, "linewidth", 3 );
    xlim ( [0.38, 0.46 ] );
    xlabel ("tGPS - tCenter [s]");
    ylabel ("Strain/SX(f)");
    title ( sprintf ( "TS-OW-%s", bname ) );
    grid on;
    hold off;

  endif

  return;
endfunction
