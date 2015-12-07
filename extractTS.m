#!/usr/bin/octave -q

function tS = extractTS ( varargin )

  uvar = parseOptions ( varargin,
                        {"fMin", "real,strictpos,scalar", 200 },
                        {"fMax", "real,strictpos,scalar", 300 },
                        {"Twindow", "real,strictpos,scalar", 10 },	%% time-window +- to extract around the event
                        {"showPlots", "bool", false },
                        {"simulate", "bool", false }
                      );
  assert ( uvar.fMax > uvar.fMin );
  tS = [];
  tEvent = 1126259462; ## from GraceDB https://gracedb.ligo.org/events/view/G184098

  fudge = 10 * eps;
  %% load frequency-domain data from SFTs:
  fNy = 1370;	%% 2x1370Hz sampling, enough to allow for resolved ~7.3ms time-shift ~ 20bins
  fnames = {"H-1_H1_1800SFT_ER8-1126257832-1800.sft"; "L-1_L1_1800SFT_ER8-1126258841-1800.sft" };

  nukefreqs = []; %% [ 120, 128, 170.5, 180, 196.1, 218.8, 240, 244.4, 256, 300, 302.2, 303, 306.2, 307.5, 315, 318.3, 331.3, 331.9, 352, 392.2, 499.6 ];	%% will generously nuke around each

  %% list of frequencies to nuke, identified "by eye" from the SFT PSD spectrum in the range [50, 500]Hz
  %% nukefreqs{1} = [ 60, 64, 73.9, 96.6, 120, 128, 144, 160, 180, 256, 288, 299.54, 299.705, 302.22, 303.305, 331.9, 352 ];	%% H1
  %% nukefreqs{2} = [ 331.3, 499.6 ]; %% L1
  nukewidth = 0.1;

  iFig = 1;
  for X = 1:length(fnames)
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
      extraLabel = "-sim";
    else
      xk0{X} = sft.SFTdata(:,1) + I * sft.SFTdata(:,2);
      extraLabel = "";
    endif
    assert ( length(fk0) == length(xk0) );

    ## %% nuke lines
    ## if ( length ( nukefreqs{X} ) > 0 )
    ##   for i = 1 : length(nukefreqs{X})
    ##     inds_nuke = find ( (fk0{X} >= nukefreqs{X}(i) - nukewidth) & (fk0{X} <= nukefreqs{X}(i) + nukewidth) );
    ##     xk0{X}(inds_nuke) = 0;
    ##   endfor
    ## endif

    %% extract frequency band of interest [fMin,fMax]
    inds0 = find ( (fk0{X} >= uvar.fMin * (1-fudge) ) & (fk0{X} <= uvar.fMax * ( 1 + fudge) ) );

    Sn{X} = median ( abs(xk0{X}(inds0)).^2 );
    xkNorm{X} = xk0{X}(inds0) ./ sqrt( Sn{X} );
    autolines_inds = find ( abs ( xkNorm{X} ) > 6 );
    autolines{X} = fk0{X} ( inds0 ( autolines_inds ) );

    ## %% nuke lines auto-identified as frequencies exceeding norm > 5
    for i = 1 : length(autolines{X})
      inds_nuke = find ( (fk0{X} >= autolines{X}(i) - nukewidth) & (fk0{X} <= autolines{X}(i) + nukewidth) );
      xk0{X}(inds_nuke) = 0;
    endfor

    if ( uvar.showPlots )
      sleg = sprintf (";%s;", IFO{X} );

      figure(iFig++); clf;
      plot ( fk0{X}(inds0), abs(xkNorm{X}), sleg );

      figure(iFig++); clf;
      plot ( fk0{X}(inds0), abs(xk0{X}(inds0)), sleg );
    endif

    %% place this band into a full spectrum including negative frequencies
    fk1 = -fNy : df : fNy;
    xk1 = zeros ( size(fk1) );
    inds1P = find ( (fk1 >= uvar.fMin * ( 1-fudge) ) & (fk1 <= uvar.fMax * ( 1 + fudge) ) );
    assert ( length(inds0) == length(inds1P) );
    xk1(inds1P) = xk0{X}(inds0);

    inds1N = find ( (fk1 >= -uvar.fMax * ( 1 + fudge) ) & (fk1 <= -uvar.fMin * ( 1 - fudge) ) );
    assert ( length(inds0) == length(inds1N) );
    xk1(inds1N) = conj ( xk0{X}( flipdim (inds0) ) );	%% mirror-image and complex-conjugate

    tS0 = FourierTransformInv ( fk1, xk1 );
    tS0.ti += t0;	%% label times by correct epoch

    inds = find ( (tS0.ti >= (tEvent - uvar.Twindow)) & (tS0.ti <= (tEvent + uvar.Twindow )) );
    tS{X}.ti = tS0.ti ( inds );
    tS{X}.xi = tS0.xi ( inds );

    %% check that we constructed a real-valued timeseries
    err = max ( abs(imag( tS{X}.xi)) ./ abs(real(tS{X}.xi) ) );
    assert ( err < 1e-6 );
    tS{X}.xi = real ( tS{X}.xi );

    out_fname = sprintf ( "TS-%s-freq%.0fHz-%.0fHz-%.0fs%s.dat", IFO{X}, uvar.fMin, uvar.fMax, 2*uvar.Twindow, extraLabel);
    fid = fopen ( out_fname, "wb" );
    fprintf ( fid, "%.9f  %g\n", [tS{X}.ti', tS{X}.xi']' );
    fclose(fid);
  endfor %% X

  if ( uvar.showPlots )
    figure(iFig++); clf; hold on;
    sleg1 = sprintf (";%s;", IFO{1} );
    sleg2 = sprintf (";%s;", IFO{2} );
    plot ( tS{1}.ti - tEvent, tS{1}.xi, sleg1, "linewidth", 3, tS{2}.ti - tEvent + 7.3e-3, (-1)*tS{2}.xi, sleg2, "linewidth", 3 );
    xlim ( [0.38, 0.46 ] );
    xlabel ("tGPS - tEvent [s]");
    ylabel ("Strain");
    tlabel = sprintf ( "TS-freq%.0fHz-%.0fHz%s", uvar.fMin, uvar.fMax, extraLabel );
    title ( tlabel );
    grid on;
    yrange = [-1.2, 1.2] * 1e-21;
    markers = [0.425, 0.427, 0.43, 0.431, 0.432, 0.433 ];
    for l = 1:length(markers)
      line ( [ 1, 1 ] * markers(l), yrange, "linestyle", "--" );
    endfor
    ylim ( yrange );
    hold off;
    plot2pdf ( tlabel );
  endif


  return;
endfunction
