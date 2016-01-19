#!/usr/bin/octave -q

global debugLevel = 1;

%% --------------------------------------------------
%% run a ringdown search as function of start-time 't0' on GW150914
%% --------------------------------------------------

%% ========== driver parameters ==========
SFTs = {"./Data/H-1_H1_1800SFT_ER8-C01-1126257832-1800.sft"; "./Data/L-1_L1_1800SFT_ER8-C01-1126258841-1800.sft" };
fSamp = 4000;	%% full sampling frequency of fmax=2kHz SFT, and conveniently such that 7.0ms timeshift between IFOs
                %% can be represented by exactly by an integer bin-shift: 7e-3 s * 4e3 Hz =  28.0 bins

plotSummary = false;
plotSpectra = false;
useTSBuffer = true;
plotPosteriors = true;

if ( !exist ("searchType" ) )
  searchType = "verify";
endif
if ( !exist ( "extraLabel" ) )
  extraLabel = "";
endif

%% ---------- Prior range defaults ----------
%%%data_FreqRange  = [ 100, 300 ]; %% avoid nasty noise stuff > 300Hz in L1
sideband = 15;	%% from running-median window 300bins * 1/(2*T)
%%data_FreqRange  = [ 15 + sideband, 2000 - sideband - 1 ];
data_FreqRange  = [ 100, 300 ];
prior_f0Range   = [ 210, 270 ];
prior_tauRange  = [ 1e-3, 20e-3 ];
prior_H         = 4e-22;	%% allow going up from 1e-22 to ~1e-21, fairly "flat" in that region
step_f0         = 0.5;
step_tau        = 0.5e-3;

switch ( searchType )
  case "verify"
    %% ---------- test-case to compare different code-versions on ----------
    tCenter = 1126259462;
    if ( !exist ( "tOffs" ) )
      tOffs = 0.43;
    endif
    tOffsStart = tOffs;
    dtOffs     = 0.0005;
    tOffsEnd   = tOffs;

    plotSpectra = true;
    useTSBuffer = false;
    plotBSGHist = false;

  case "onSource"
    %% ---------- "ON-SOURCE ----------
    tCenter = 1126259462;
    tOffsStart = 0.422;		%% this is about when the signal enters f0>~200Hz
    dtOffs     = 0.0005;	%% 0.5ms stepsize
    tOffsEnd   = 0.435; %% 0.438;

    plotBSGHist = false;
    plotSummary = true;
    plotPosteriors = true;
    useTSBuffer = true;

  case "offSource"
    %% ---------- "OFF-SOURCE" for background estimation ----------
    tCenter     = 1126259472; %% 10s past GW150914
    tOffsStart  = -5;
    dtOffs      = 0.1;	%% 100ms stepsize, to avoid template overlap --> 'independent' templates
    tOffsEnd    = 5;
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
  [ts{X}, ft{X}, psd{X}] = extractTSfromSFT ( "SFTpath", SFTs{X}, "fMin", min(data_FreqRange), "fMax", max(data_FreqRange), "fSamp", fSamp, "tCenter", tCenter, "plotSpectrum", plotSpectra, "useBuffer", useTSBuffer );
endfor

%% create unique time-tagged 'ResultsDir' for each run:
gm = gmtime ( time () );
resDir = sprintf ( "Results/Results-%02d%02d%02d-%02dh%02d-%s-data%.0fHz-%.0fHz%s", gm.year - 100, gm.mon + 1, gm.mday, gm.hour, gm.min, searchType, data_FreqRange, extraLabel );
[status, msg, id] = mkdir ( resDir ); assert ( status == 1, "Failed to created results dir '%s': %s\n", resDir, msg );
addpath ( pwd() );
cd ( resDir );

try
%% ----- run search
tOffsV = [ tOffsStart : dtOffs : tOffsEnd ];
Nsteps = length(tOffs);

ret = cell ( 1, Nsteps );
plot = [];
for i = 1:Nsteps
  DebugPrintf ( 1, "tOffs = %.4f s:\n", tOffsV(i) );
  ret{i} = searchRingdown ( "ts", ts, "psd", psd, "tOffs", tOffsV(i), "tCenter", tCenter, ...
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
  figure(); clf;
  xrange = [ tOffsStart, tOffsEnd ];

  subplot ( 2, 2, 1, "align" );
  plot ( tOffsV, log10(dat.BSG_mean), "-o" ); grid on;
  xlim ( xrange );
  ylim ( [ -2, 5 ] );
  line ( xrange, 0, "linestyle", "-", "linewidth", 3 );
  line ( xrange, 1, "linestyle", ":", "linewidth", 3 );
  xlim ( xrange );
  ylabel ("log10<BSG>");

  subplot ( 2, 2, 3, "align" );
  plot ( tOffsV, dat.SNR_MPE, "-o" ); grid on;
  xlim ( xrange );
  ylabel ("SNR(MPE)");

  subplot ( 2, 2, 2, "align" );
  errorbar ( tOffsV, dat.f0_MPE, dat.f0_lerr, dat.f0_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( prior_f0Range );
  ylabel ("f0 [Hz]");

  subplot ( 2, 2, 4, "align" );
  errorbar ( tOffsV, 1e3*dat.tau_MPE, 1e3*dat.tau_lerr, 1e3*dat.tau_uerr, ";90%;" ); grid on;
  xlim ( xrange );
  ylim ( 1e3 * prior_tauRange );
  ylabel ("tau [ms]");
  xlabel ( sprintf ( "%.0f + tOffs [s]", tCenter) );
  fname = sprintf ( "%s-summary.pdf", ret{1}.bname );
  ezprint ( fname, "width", 512 );
endif

if ( plotBSGHist )
  figure(); clf;
  hist ( dat.BSG_mean, 20 );
  xlabel ( "<BSG>" );
  fname = sprintf ( "%s-hist.pdf", ret{1}.bname );
  ezprint ( fname, "width", 512 );
endif

cd ("../..");

catch
  err = lasterror();
  warning(err.identifier, err.message);
  err.stack.file
  err.stack.name
  err.stack.line

  cd ("../..");	%% make sure we end up in main dir in case of failure
end_try_catch
