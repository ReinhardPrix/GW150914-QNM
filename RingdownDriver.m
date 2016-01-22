#!/usr/bin/octave -q

global debugLevel = 1;

%% --------------------------------------------------
%% run a ringdown search as function of start-time 't0' on GW150914
%% --------------------------------------------------

%% ========== driver parameters ==========
SFTs = {"./Data/H-1_H1_1800SFT_ER8-C01-1126257832-1800.sft"; "./Data/L-1_L1_1800SFT_ER8-C01-1126258841-1800.sft" };
fSamp = 4000;	%% full sampling frequency of fmax=2kHz SFT, and conveniently such that 7.0ms timeshift between IFOs
                %% can be represented by exactly by an integer bin-shift: 7e-3 s * 4e3 Hz =  28.0 bins
tEvent = 1126259462;
tMergerOffs = 0.42285; %% https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TestingGR/O1/G184098/ringdown_presence

plotSummary = false;
plotSpectra = false;
useTSBuffer = true;
plotPosteriors = true;

if ( !exist("searchType") )     searchType = "verify"; endif
if ( !exist("extraLabel") )     extraLabel = ""; endif
if ( !exist("psd_version") )    global psd_version = 2; endif
if ( !exist("cleanLines") )     global cleanLines = false; endif
if ( !exist("data_FreqRange") ) data_FreqRange  = [ 30, 1e3 ]; endif
if ( !exist("iFig0") )          global iFig0 = 0; endif

%% ---------- Prior range defaults ----------
prior_f0Range   = [ 200, 300 ];
prior_tauRange  = [ 0.5e-3, 20e-3 ];
prior_H         = 4e-22;	%% allow going up from 1e-22 to ~1e-21, fairly "flat" in that region
step_f0         = 0.5;
step_tau        = 0.5e-3;

%% NR best-fit 'merger time'
%% from https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/TestingGR/O1/G184098/ringdown_presence

switch ( searchType )
  case "verify"
    %% ---------- test-case to compare different code-versions on ----------
    tCenter = tEvent;
    if ( !exist ( "tOffs" ) )
      tOffs = 0.42985;
    endif
    tOffsStart = tOffs;
    dtOffs     = 0.0005;
    tOffsEnd   = tOffs;

    plotSpectra = true;
    useTSBuffer = false;
    plotBSGHist = false;

  case "onSource"
    %% ---------- "ON-SOURCE ----------
    tCenter = tEvent;
    tOffsStart = tMergerOffs;		%% this is about when the signal enters f0>~200Hz
    dtOffs     = 0.0005;	%% 0.5ms stepsize
    tOffsEnd   = tMergerOffs + 10e-3;

    plotBSGHist = false;
    plotSummary = true;
    plotPosteriors = true;
    useTSBuffer = true;

  case "offSource"
    %% ---------- "OFF-SOURCE" for background estimation ----------
    tCenter     = tEvent + 10;
    tOffsStart  = -3;
    dtOffs      = 0.05;	%% stepsize to avoid template overlap --> 'independent' templates
    tOffsEnd    = 3;
    plotPosteriors = true;
    plotSummary = true;
    plotSpectra = false;
    useTSBuffer = true;
    plotBSGHist = true;
  otherwise
    error ("Unknown searchType = '%s' specified\n", searchType );

endswitch

%% ----- data preparation -----

%% load frequency-domain data from SFTs:
for X = 1:length(SFTs)
  [ts{X}, ft{X}, psd{X}] = extractTSfromSFT ( "SFTpath", SFTs{X}, "fMin", min(data_FreqRange), "fMax", max(data_FreqRange), "fSamp", fSamp, ...
                                              "tCenter", tCenter, "Twindow", 4, ...
                                              "plotSpectrum", plotSpectra, "useBuffer", useTSBuffer );
endfor

%% create unique time-tagged 'ResultsDir' for each run:
gm = gmtime ( time () );
resDir = sprintf ( "Results/Results-%02d%02d%02d-%02dh%02d-%s-data%.0fHz-%.0fHz-psd_v%d-lineCleaning%s%s",
                   gm.year - 100, gm.mon + 1, gm.mday, gm.hour, gm.min, searchType, data_FreqRange,
                   psd_version, ifelse ( cleanLines, "On", "Off" ),
                   extraLabel );
[status, msg, id] = mkdir ( resDir ); assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );
addpath ( pwd() );
cd ( resDir );

%%try
%% ----- run search
tOffsV = [ tOffsStart : dtOffs : tOffsEnd ];
Nsteps = length(tOffsV);

ret = cell ( 1, Nsteps );
for i = 1:Nsteps
  DebugPrintf ( 1, "tOffs = %.5f s:\n", tOffsV(i) );
  ret{i} = searchRingdown ( "ts", ts, "psd", psd, ...
                            "tCenter", tCenter, "tOffs", tOffsV(i), ...
                            "prior_f0Range", prior_f0Range, "step_f0", step_f0, ...
                            "prior_tauRange", prior_tauRange, "step_tau", step_tau, ...
                            "prior_H", prior_H, ...
                            "plotResults", plotPosteriors
                          );

  %% ----- save posterior in matrix format ----------
  if ( i == 1 )
    fname = sprintf ( "f0s.dat", ret{i}.bname );
    tmp = ret{i}.ff0;
    save ( "-ascii", fname, "tmp" );
    fname = sprintf ( "taus.dat", ret{i}.bname );
    tmp = ret{i}.ttau;
    save ( "-ascii", fname, "tmp" );
  endif
  fname = sprintf ( "%s-BSG.dat", ret{i}.bname );
  tmp = ret{i}.BSG;
  save ( "-ascii", fname, "tmp" );

  fname = sprintf ( "%s-SNR.dat", ret{i}.bname );
  tmp = ret{i}.SNR;
  save ( "-ascii", fname, "tmp" );

  %% ---------- collect results for plotting
  dat.BSG_mean(i) = ret{i}.BSG_mean;
  dat.A_MPE(i)    = ret{i}.A_MPE;
  dat.phi0_MPE(i) = ret{i}.phi0_MPE;
  dat.SNR_MPE(i)  = ret{i}.SNR_MPE;

  dat.f0_MPE(i)   = ret{i}.f0_est.MPE;
  dat.f0_lerr(i)  = dat.f0_MPE(i) - ret{i}.f0_est.lower;
  dat.f0_uerr(i)  = ret{i}.f0_est.upper - dat.f0_MPE(i);

  dat.tau_MPE(i)   = ret{i}.tau_est.MPE;
  dat.tau_lerr(i)  = dat.tau_MPE(i) - ret{i}.tau_est.lower;
  dat.tau_uerr(i)  = ret{i}.tau_est.upper - dat.tau_MPE(i);

  %% ---------- plot results summary page
  %%fname = sprintf ( "%s.png", ret{i}.bname );
  %%ezprint ( fname, "width", 1024, "height", 786, "dpi", 72 );
  if ( plotPosteriors )
    fname = sprintf ( "%s.pdf", ret{i}.bname );
    ezprint ( fname, "width", 512 );
  endif

endfor

%% ----- store complete results dump ----------
fname = sprintf ( "RingdownDriver-%s.hd5", ret{1}.bname );
save ("-hdf5", fname )

%% ----- plot quantities vs tOffs ----------
if ( plotSummary )
  figure ( iFig0 * 5 + 4 ); clf;

  %% ----- plot log10BSG(tOffs)
  subplot ( 2, 2, 1, "align" );
  xrange = [ tOffsStart, tOffsEnd ];
  yrange = [ -1, 10 ];
  plot ( tOffsV, log10(dat.BSG_mean), "-o" );
  xlim ( xrange );
  ylim ( yrange );
  line ( xrange, 0, "linestyle", "-", "linewidth", 3 );
  line ( xrange, 1, "linestyle", ":", "linewidth", 3 );
  line ( tMergerOffs * [1,1], yrange, "linestyle", "-", "linewidth", 2 );
  xlim ( xrange );
  ylabel ("log10<BSG>");
  xlabel ( "tOffs [s]")
  grid on;
  %% add second x-axis on top
  axes1 = gca ();
  set (axes1, "XAxisLocation",  "bottom");
  set (axes1, "activepositionproperty", "position")
  hold on;
  axes2 = axes ();
  set (axes2, "color", "none", "ytick", [])
  set (axes2, "XAxisLocation",  "top" )
  set (axes2, "activepositionproperty", "position")
  set (axes2, "position", get (axes1, "position"))
  set (axes2, "XTick", (get ( axes1, "xtick" ) - tMergerOffs) * 1e3 );
  hold off
  set (axes2, "xlim", (get (axes1, "xlim") - tMergerOffs) * 1e3 );
  xlabel ( "tOffs from Merger [ms]")

  %% ----- plot f0_MPE(tOffs)
  subplot ( 2, 2, 2, "align" );
  errorbar ( tOffsV, dat.f0_MPE, dat.f0_lerr, dat.f0_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( prior_f0Range );
  ylabel ("f0 [Hz]");
  %% add second x-axis on top
  axes1 = gca ();
  set (axes1, "XAxisLocation",  "bottom");
  set (axes1, "activepositionproperty", "position")
  hold on;
  axes2 = axes ();
  set (axes2, "color", "none", "ytick", [])
  set (axes2, "XAxisLocation",  "top" )
  set (axes2, "activepositionproperty", "position")
  set (axes2, "position", get (axes1, "position"))
  set (axes2, "XTick", (get ( axes1, "xtick" ) - tMergerOffs) * 1e3 );
  hold off
  set (axes2, "xlim", (get (axes1, "xlim") - tMergerOffs) * 1e3 );
  xlabel ( "tOffs from Merger [ms]")


  %% ----- plot SNR(tOffs)
  subplot ( 2, 2, 3, "align" );
  plot ( tOffsV, dat.SNR_MPE, "-o" ); grid on;
  xlim ( xrange );
  ylabel ("SNR(MPE)");
  %% add second x-axis on top
  axes1 = gca ();
  set (axes1, "XAxisLocation",  "bottom");
  set (axes1, "activepositionproperty", "position")
  hold on;
  axes2 = axes ();
  set (axes2, "color", "none", "ytick", [])
  set (axes2, "XAxisLocation",  "top" )
  set (axes2, "activepositionproperty", "position")
  set (axes2, "position", get (axes1, "position"))
  set (axes2, "XTick", (get ( axes1, "xtick" ) - tMergerOffs) * 1e3 );
  hold off
  set (axes2, "xlim", (get (axes1, "xlim") - tMergerOffs) * 1e3 );


  %% ----- plot tau_MPE(tOffs)
  subplot ( 2, 2, 4, "align" );
  errorbar ( tOffsV, 1e3*dat.tau_MPE, 1e3*dat.tau_lerr, 1e3*dat.tau_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( 1e3 * prior_tauRange );
  ylabel ("tau [ms]");
  xlabel ( sprintf ( "%.0f + tOffs [s]", tCenter) );
  %% add second x-axis on top
  axes1 = gca ();
  set (axes1, "XAxisLocation",  "bottom");
  set (axes1, "activepositionproperty", "position")
  hold on;
  axes2 = axes ();
  set (axes2, "color", "none", "ytick", [])
  set (axes2, "XAxisLocation",  "top" )
  set (axes2, "activepositionproperty", "position")
  set (axes2, "position", get (axes1, "position"))
  set (axes2, "XTick", (get ( axes1, "xtick" ) - tMergerOffs) * 1e3 );
  hold off
  set (axes2, "xlim", (get (axes1, "xlim") - tMergerOffs) * 1e3 );


  fname = sprintf ( "%s-summary.pdf", ret{1}.bname );
  ezprint ( fname, "width", 512 );

endif

if ( plotBSGHist )
  figure ( iFig0 * 5  + 5 ); clf;
  hist ( dat.BSG_mean, 20 );
  xlabel ( "<BSG>" );
  fname = sprintf ( "%s-hist.pdf", ret{1}.bname );
  ezprint ( fname, "width", 512 );
endif

cd ("../..");

## catch
##   err = lasterror();
##   warning(err.identifier, err.message);
##   err.stack.file
##   err.stack.name
##   err.stack.line

##   cd ("../..");	%% make sure we end up in main dir in case of failure
## end_try_catch
